package multicell

import (
	"gonum.org/v1/gonum/mat"
)

func NewTensor3(n0, n1, n2 int) Tensor3 {
	ten := make([]Dmat, n0)
	for i := range ten {
		mat := make([]Vec, n1)
		for i := range mat {
			mat[i] = NewVec(n2)
		}
		ten[i] = mat
	}
	return ten
}

func flatten0(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len0)

	for i := 0; i < len0; i++ {
		for j := 0; j < len1; j++ {
			for k := 0; k < len2; k++ {
				mat[i] = append(mat[i], ten[i][j][k])
			}
		}
	}

	return mat
}

func flatten1(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len1)

	for j := 0; j < len1; j++ {
		for i := 0; i < len0; i++ {
			for k := 0; k < len2; k++ {
				mat[j] = append(mat[j], ten[i][j][k])
			}
		}
	}

	return mat
}

func flatten2(ten Tensor3) Dmat {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat := make([]Vec, len2)

	for k := 0; k < len2; k++ {
		for i := 0; i < len0; i++ {
			for j := 0; j < len1; j++ {
				mat[k] = append(mat[k], ten[i][j][k])
			}
		}
	}

	return mat
}

func FlattenTensor3(ten Tensor3, ind int) Dmat {
	if ind == 0 {
		return flatten0(ten)
	} else if ind == 1 {
		return flatten1(ten)
	} else {
		return flatten2(ten)
	}
}

func GetHOSVD(ten Tensor3) (Tensor3, *mat.Dense, *mat.Dense, *mat.Dense) {
	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat0 := flatten0(ten)
	u0, s0, _ := GetSVD(mat0)
	mat1 := flatten1(ten)
	u1, s1, _ := GetSVD(mat1)
	mat2 := flatten2(ten)
	u2, s2, _ := GetSVD(mat2)

	rank0 := len(s0)
	rank1 := len(s1)
	rank2 := len(s2)

	score := NewTensor3(len0, len1, len2)
	for i := 0; i < len0; i++ {
		for j := 0; j < len1; j++ {
			for k := 0; k < len2; k++ {
				go func(i, j, k int) {
					for i1 := 0; i1 < rank0; i1++ {
						for j1 := 0; j1 < rank1; j1++ {
							for k1 := 0; k1 < rank2; k1++ {
								score[i][j][k] += ten[i][j][k] * u0.At(i, i1) * u1.At(j, j1) * u2.At(k, k1)
							}
						}
					}
				}(i, j, k)
			}
		}
	}
	return score, u0, u1, u2
}

// Interlaced Computation (Sequentially Truncated HOSVD); supposed to be significantly faster.
func GetST_HOSVD(ten Tensor3) (Tensor3, *mat.Dense, *mat.Dense, *mat.Dense) {

	len0 := len(ten)
	len1 := len(ten[0])
	len2 := len(ten[0][0])

	mat0 := flatten0(ten)
	u0, s0, _ := GetSVD(mat0)
	rank0 := len(s0)
	A0 := NewTensor3(len0, len1, len2)
	for n := 0; n < rank0; n++ {
		u := u0.ColView(n)
		for i := 0; i < len0; i++ {
			ui := u.AtVec(i)
			for j := 0; j < len1; j++ {
				for k := 0; k < len2; k++ {
					A0[i][j][k] += ten[i][j][k] * ui
				}
			}
		}
	}

	mat1 := flatten1(A0)
	u1, s1, _ := GetSVD(mat1)
	rank1 := len(s1)
	A1 := NewTensor3(len0, len1, len2)
	for n := 0; n < rank1; n++ {
		u := u1.ColView(n)
		for j := 0; j < len1; j++ {
			uj := u.AtVec(j)
			for i := 0; i < len0; i++ {
				for k := 0; k < len2; k++ {
					A1[i][j][k] += A0[i][j][k] * uj
				}
			}
		}
	}

	mat2 := flatten2(A1)
	u2, s2, _ := GetSVD(mat2)
	rank2 := len(s2)
	A2 := NewTensor3(len0, len1, len2)
	for n := 0; n < rank2; n++ {
		u := u2.ColView(n)
		for k := 0; k < len1; k++ {
			uk := u.AtVec(k)
			for i := 0; i < len0; i++ {
				for j := 0; j < len1; j++ {
					A2[i][j][k] += A1[i][j][k] * uk
				}
			}
		}
	}

	return A2, u0, u1, u2
}

func GetCrossCov3(vecs0, vecs1, vecs2 []Vec, submean0, submean1, submean2 bool) (Vec, Vec, Vec, Tensor3) {
	cv0 := GetMeanVec(vecs0)
	cv1 := GetMeanVec(vecs1)
	cv2 := GetMeanVec(vecs2)
	len0, len1, len2 := len(cv0), len(cv1), len(cv2)

	cov := NewTensor3(len0, len1, len2)

	var d0, d1, d2 float64

	for n := range vecs0 {
		for i, c0 := range cv0 {
			if submean0 {
				d0 = vecs0[n][i] - c0
			} else {
				d0 = vecs0[n][i]
			}
			for j, c1 := range cv1 {
				if submean1 {
					d1 = vecs1[n][j] - c1
				} else {
					d1 = vecs1[n][j]
				}
				for k, c2 := range cv2 {
					if submean2 {
						d2 = vecs2[n][k] - c2
					} else {
						d2 = vecs2[n][k]
					}
					cov[i][j][k] += d0 * d1 * d2
				}
			}
		}
	}

	fn := 1.0 / float64(len(vecs0))
	for i := range cv0 {
		for j := range cv1 {
			for k := range cv2 {
				cov[i][j][k] *= fn
			}
		}
	}

	return cv0, cv1, cv2, cov
}
