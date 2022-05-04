package multicell

// Canonical Polyadic Decomposition of 3-tensors.
/*
   The ALS1-CPO Algorithm in
   Mikael Sorensen, Lieven de Lathauwer, Pierre Comon, Sylvie Icart, Luc Deneire.
   Canonical Polyadic decomposition with a Columnwise Orthonormal Factor Matrix.
   SIAM Journal on Matrix Analysis and Applications, 2012, 33 (4), pp.1190-1213.
   ( https:/doi.org/10.1137/110830034 )
*/

import (
	"fmt"
	"log"
	"math"
	//	"gonum.org/v1/gonum/mat"
)

// X_{(1)}^T
func tflatten0(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len2) // (len0*len1) X len2
	for k := 0; k < len2; k++ {
		for i := 0; i < len0; i++ {
			for j := 0; j < len1; j++ {
				mat[k] = append(mat[k], ten[i][j][k])
			}
		}
	}
	return mat
}

func tflatten1(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len0) // (len1*len2) X len0

	for i := 0; i < len0; i++ {
		for j := 0; j < len1; j++ {
			for k := 0; k < len2; k++ {
				mat[i] = append(mat[i], ten[i][j][k])
			}
		}
	}
	return mat
}

func tflatten2(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len1) // (len2*len0) X len1

	for j := 0; j < len1; j++ {
		for k := 0; k < len2; k++ {
			for i := 0; i < len0; i++ {
				mat[j] = append(mat[j], ten[i][j][k])
			}
		}
	}
	return mat
}

// The Khatri-Rao Product between KxN and KxN matrices
func KR_Product(a, b Dmat) Dmat {
	lena := len(a)
	lenb := len(b)
	lenc := len(a[0]) // should be equal to len(b[0])

	mat := NewDmat(lena*lenb, lenc)

	for i := 0; i < lena; i++ {
		ai := a[i]
		for j := 0; j < lenb; j++ {
			ind := i*lenb + j
			bj := b[j]
			for k := 0; k < lenc; k++ {
				mat[ind][k] = ai[k] * bj[k]
			}
		}
	}

	return mat
}

func NormDiag2(a Dmat) Vec {
	len0 := len(a)
	rank := len(a[0])
	dvec := NewVec(rank)

	for j := 0; j < rank; j++ {
		d := 0.0
		for i := 0; i < len0; i++ {
			v := a[i][j]
			d += v * v
		}
		if d > 0.0 {
			dvec[j] = 1.0 / d
		} else {
			dvec[j] = 0.0
		}
	}

	return dvec
}

// Canonical Polyadic Decomposition
// X = S*(A1 x A2 x A3) where A3 is orthogonal.
func GetCPDO(ten Tensor3) (Vec, Dmat, Dmat, Dmat) {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	rank := len0
	if len1 < rank {
		rank = len1
	}
	if len2 < rank {
		rank = len2
	}

	a0 := NewDmat(len0, rank)
	a1 := NewDmat(len1, rank)
	a2 := NewDmat(len2, rank)
	ta2 := NewDmat(len2, rank)

	// initialize to the "Identity" matrix.
	for i := 0; i < rank; i++ {
		a0[i][i] = 1.0
		a1[i][i] = 1.0
		a2[i][i] = 1.0
	}
	xt0 := tflatten0(ten) // len2 x (len0*len1)
	xt1 := tflatten1(ten) // len0 x (len1*len2)
	tx2 := tflatten2(ten) // len1 x (len2*len0)
	dev := 1000.0

	for istep := 0; istep < 100; istep++ {
		pa01 := KR_Product(a0, a1)   // (len0*len1) x rank
		tmat := MultDmats(xt0, pa01) // len2 x rank

		// Step 1: Update a0.
		U, sval, V := GetSVD(tmat) // U: len2 x rank; V: rank x rank
		log.Println("GetCPDO Step 1", istep)

		// Step 2: update a2.
		ResetDmat(ta2)
		for i := 0; i < len2; i++ {
			for j := range sval {
				for k := range sval {
					ta2[i][k] += U.At(i, j) * V.At(k, j)
				}
			}
		}
		log.Println("GetCPDO Step 2", istep)
		// check convergence
		dev = 0.0
		mag := 0.0
		for i := 0; i < len2; i++ {
			for j := range sval {
				d := ta2[i][j] - a2[i][j]
				mag += ta2[i][j] * ta2[i][j]
				dev += d * d
			}
		}
		dev /= mag
		CopyDmat(a2, ta2)

		// Step 3: Update a0.
		pa12 := KR_Product(a1, a2)  // (len1*len2) x rank
		ta0 := MultDmats(xt1, pa12) // len0 x rank
		d1 := NormDiag2(a1)         // rank
		for i := 0; i < len0; i++ {
			for j := 0; j < rank; j++ {
				a0[i][j] = ta0[i][j] * d1[j]
			}
		}
		log.Println("GetCPDO Step 3", istep)

		// Step 4: Normalize a0.
		d0 := NormDiag2(a0)
		for i := 0; i < len0; i++ {
			for j := 0; j < rank; j++ {
				a0[i][j] *= math.Sqrt(d0[j])
			}
		}
		log.Println("GetCPDO Step 4", istep)

		// Step 5: Update a1.
		pa20 := KR_Product(a2, a0)  // (len2*len0) x rank
		ta1 := MultDmats(tx2, pa20) // len1 x rank
		CopyDmat(a1, ta1)

		log.Println("GetCPDO Step 5", istep)
		log.Println("GetCPDO:", istep, dev)
		fmt.Println("A0", a0)
		fmt.Println("A1", a1)
		fmt.Println("A2", a2)

	}

	sigma := NewVec(rank)
	for _, ai := range a1 {
		for j, aij := range ai {
			sigma[j] += aij * aij
		}
	}

	for i, s := range sigma {
		sigma[i] = math.Sqrt(s)
	}

	return sigma, a0, a1, a2
}
