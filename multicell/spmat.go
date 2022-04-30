package multicell

import (
	//	"math"
	"math/rand"
)

type Spmat struct {
	Ncol int                 // number of columns
	Mat  [](map[int]float64) // Sparse matrix is an array of maps.
}

func NewSpmat(nrow, ncol int) Spmat { //Initialize new sparse matrix
	mat := make([](map[int]float64), nrow)
	for i := range mat {
		mat[i] = make(map[int]float64)
	}
	return Spmat{ncol, mat}
}

func (sp *Spmat) Copy() Spmat {
	nsp := NewSpmat(len(sp.Mat), sp.Ncol)
	for i, m := range sp.Mat {
		for j, d := range m {
			nsp.Mat[i][j] = d
		}
	}
	return nsp
}

func (sp *Spmat) Randomize(density, sd float64) { //Randomize entries of sparse matrix
	for i := range sp.Mat {
		for j := 0; j < sp.Ncol; j++ {
			r := rand.Float64()
			if r < density {
				sp.Mat[i][j] = rand.NormFloat64() * sd //Scale to theoretical sd per entry
			}
		}
	}
}

func DiffSpmat(m1, m2 *Spmat) Spmat { //This function works fine
	d := NewSpmat(len(m1.Mat), m1.Ncol) //initialization
	ncol := m1.Ncol
	for i := range m1.Mat {
		for j := 0; j < ncol; j++ {
			d.Mat[i][j] = m1.Mat[i][j] - m2.Mat[i][j]
		}
	}

	return d
}

func (sp *Spmat) Scale(c float64) {
	for i, mi := range sp.Mat {
		for j := range mi {
			sp.Mat[i][j] *= c
		}
	}
}

func MultMatVec(vout Vec, mat Spmat, vin Vec) { //Matrix multiplication
	for i := range vout {
		vout[i] = 0.0
	}

	for i, m := range mat.Mat {
		for j, d := range m {
			vout[i] += d * vin[j]
		}
	}
	return
}

func MultMatVec_T(vout Vec, mat Spmat, vin Vec) { //Matrix transposition and then multiplication
	for i := range vout {
		vout[i] = 0.0
	}
	for i, m := range mat.Mat {
		vi := vin[i]
		for j, d := range m {
			vout[j] += d * vi
		}
	}

	return
}

func (mat *Spmat) mutateSpmat(density, sd float64) { //mutating a sparse matrix
	nrow := len(mat.Mat)
	nmut := int(mutRate * float64(nrow*mat.Ncol))
	for n := 0; n < nmut; n++ {
		i := rand.Intn(nrow)
		j := rand.Intn(mat.Ncol)
		r := rand.Float64()
		delete(mat.Mat[i], j)
		if r < density {
			mat.Mat[i][j] = rand.NormFloat64() * sd //Scale to theoretical sd per entry.
		}
	}
	//Note: This implementation has non-zero probability of choosing same element to be mutated twice.
	return
}

func DotSpmats(mat0, mat1 Spmat) float64 {
	dot := 0.0
	for i, m := range mat0.Mat {
		for j, d := range m {
			dot += d * mat1.Mat[i][j]
		}
	}

	return dot
}
