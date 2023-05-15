package multicell

import (
	"fmt"
	"log"
	"math"
	"math/rand"

	"gonum.org/v1/gonum/mat"
)

var maxPop int = 200 // population size
var maxDevStep = 200 // Maximum steps for development.
var ngenes int = 200 // number of genes
var nenv int = 200   // number of environmental cues/traits per face
var nsel int = 40    // number of environmental cues/traits per face FOR SELECTION
var ncells int = 1   //number of cell types

var with_cue bool = true // with or without environmental cues.
var pheno_feedback bool = false
var withE bool = false // = with_cue || pheno_feedback
var withF bool = true  // Epigenetic marker layer
var withH bool = true  // Higher order complexes layer
var withJ bool = false
var devNoise float64 = 0.05 // noise strength
var mutRate float64 = 0.005 // mutation rate

// Decay rates
var tauF float64 = 0.2
var tauG float64 = 1.0
var tauH float64 = 1.0

// Damping rate of environmental cues
var dampFactorE float64 = 1.0

// matrix densities
const defaultDensity float64 = 0.02 // = 4 / 200 (4 inputs per row)
var DensityE float64 = defaultDensity
var DensityF float64 = defaultDensity
var DensityG float64 = defaultDensity
var DensityH float64 = defaultDensity
var DensityJ float64 = defaultDensity
var DensityP float64 = defaultDensity

type Settings struct {
	MaxPop     int // Maximum number of individuals in population
	MaxDevStep int // Maximum steps for development.
	NGenes     int
	NEnv       int
	NSel       int
	NCells     int     // Number of cell types
	WithCue    bool    // With cue?
	FLayer     bool    // f present?
	HLayer     bool    // h present?
	JLayer     bool    //  J present?
	Pfback     bool    // P feedback to E layer
	SDNoise    float64 // probability (or stdev) of environmental noise
	MutRate    float64 // mutation rate
	TauF       float64
	TauG       float64
	TauH       float64
	DensityE   float64
	DensityF   float64
	DensityG   float64
	DensityH   float64
	DensityJ   float64
	DensityP   float64
}

func CurrentSettings() Settings {
	return Settings{MaxPop: maxPop, MaxDevStep: maxDevStep,
		NGenes: ngenes, NEnv: nenv, NSel: nsel, NCells: ncells,
		WithCue: with_cue, FLayer: withF, HLayer: withH, JLayer: withJ,
		Pfback: pheno_feedback, SDNoise: devNoise, MutRate: mutRate,
		TauF: tauF, TauG: tauG, TauH: tauH,
		DensityE: DensityE, DensityF: DensityF, DensityG: DensityG,
		DensityH: DensityH, DensityJ: DensityJ, DensityP: DensityP}

}

const (
	cueMag = 1.0 // each trait is +/-cueMag
	//maxDevStep = 200    // Maximum steps for development.
	epsDev   = 1.0e-5 // Convergence criterion of development.
	eps      = 1.0e-50
	sqrt3    = 1.73205080756887729352744634150587236694280525381038062805580697
	ccStep   = 5.0                  // Number of steady steps for convergence
	alphaEMA = 2.0 / (1.0 + ccStep) // exponential moving average/variance
)

// Length of a gene for Unicellular organism.
var fullGeneLength = 4*ngenes + 2*nenv

//calculated from layers present or absent.

const baseSelStrength float64 = 20.0 // default selection strength; to be normalized by number of cells
const selDevStep float64 = 20.0      // Developmental steps for selection

const minWagnerFitness float64 = 0.01

// slope of activation functions
var omega_f float64 = 1.0
var omega_g float64 = 1.0
var omega_h float64 = 1.0
var omega_p float64 = 1.0

//Remark: defaults to full model!

type Vec = []float64 //Vector is a slice
type Dmat = []Vec
type Tensor3 []Dmat

func SetSeed(seed int64) {
	rand.Seed(seed)
}

func SetParams(s Settings) {
	maxPop = s.MaxPop
	maxDevStep = s.MaxDevStep
	ngenes = s.NGenes
	nenv = s.NEnv
	nsel = s.NSel
	ncells = s.NCells
	with_cue = s.WithCue
	withF = s.FLayer
	withH = s.HLayer
	withJ = s.JLayer
	pheno_feedback = s.Pfback
	devNoise = s.SDNoise
	mutRate = s.MutRate
	withE = with_cue || pheno_feedback
	tauF = s.TauF
	tauG = s.TauG
	tauH = s.TauH
	DensityE = s.DensityE
	DensityF = s.DensityF
	DensityG = s.DensityG
	DensityH = s.DensityH
	DensityJ = s.DensityJ
	DensityP = s.DensityP

	from_g := DensityG * float64(ngenes)

	if withE {
		from_e := DensityE * float64(nenv)
		if with_cue && pheno_feedback {
			omega_f = 1.0 / math.Sqrt(2*from_e+from_g*(2-tauG))
		} else {
			omega_f = 1.0 / math.Sqrt(from_e+from_g*(2-tauG))
		}
	} else {
		omega_f = 1.0 / math.Sqrt(from_g*(2-tauG))
		DensityE = 0.0
	}

	if withF {
		omega_g = 1.0 / math.Sqrt(DensityF*float64(ngenes)*(2-tauF))
	} else {
		DensityF = 0.0
		efac := 1.0
		if with_cue && pheno_feedback {
			efac = 2.0
		}
		omega_g = 1.0 / math.Sqrt(DensityG*float64(ngenes)*(2-tauG)+efac*DensityE*float64(nenv))
	}

	if withH {
		if withJ {
			omega_h = 1.0 / math.Sqrt(from_g*((2-tauG)+(2-tauH)))
		} else {
			omega_h = 1.0 / math.Sqrt(from_g*(2-tauG))
			DensityJ = 0.0
		}
		omega_p = 1.0 / math.Sqrt(DensityP*float64(ngenes)*(2-tauH))
	} else {
		omega_h = 0.0
		DensityH = 0.0
		omega_p = 1.0 / math.Sqrt(DensityP*float64(ngenes)*(2-tauG))
	}
	/* // Trying not calling EMA for NoDev instead.
	if maxDevStep == 1 {
		omega_p *= 10.0 //Arbitrary factor to increase sensitivity of NoDev.
	}
	*/

}

func GetMaxPop() int {
	return maxPop
}

func GetNcells() int {
	return ncells
}

func GetNenv() int {
	return nenv
}

func GetNsel() int {
	return nsel
}

func sigmoid(x, omega float64) float64 {
	return 1 / (1 + math.Exp(-x/omega))
}

func tanh(x, omega float64) float64 {
	return math.Tanh(omega * x)
}

func lecuntanh(x float64) float64 { //Le'Cun's hyperbolic tangent
	return 1.7159 * math.Tanh(2*x/3)
}

func arctan(x, omega float64) float64 {
	return math.Atan(omega * x)
}

func lecunatan(x float64) float64 { //Rescaled arctan under same treatment of Le'Cun's hyperbolic tangent.
	return 6.0 * math.Atan(x/sqrt3) / math.Pi
}

func relu(x, omega float64) float64 {
	if x < 0 {
		return 0
	} else {
		return omega * x
	}
}

func sigmaf(x float64) float64 { //Activation function for epigenetic markers
	return lecunatan(x * omega_f)
	//return tanh(x, omega_f)
}

func sigmag(x float64) float64 { //Activation function for gene expression levels
	return lecunatan(x * omega_g)
	//return tanh(x, omega_g)
}

func sigmah(x float64) float64 { //Activation function for higher order complexes
	return lecunatan(x * omega_h) //abstract level of amount of higher order complexes
	//return tanh(x, omega_h)
}

func rho(x float64) float64 { //Function for converting gene expression into phenotype
	//return cueMag * lecunatan(x*omega_p)
	return cueMag * tanh(x, omega_p)
}

func NewDmat(nrow, ncol int) Dmat {
	mat := make([]Vec, nrow)
	for i := range mat {
		mat[i] = NewVec(ncol)
	}

	return mat
}

func CopyDmat(mat1, mat0 Dmat) {
	for i, di := range mat0 {
		for j, dij := range di {
			mat1[i][j] = dij
		}
	}
}

func ResetDmat(mat Dmat) {
	for i, mi := range mat {
		for j := range mi {
			mat[i][j] = 0
		}
	}
}

func MultDmats(mat0, mat1 Dmat) Dmat {
	dimi := len(mat0)
	dimj := len(mat0[0])
	dimk := len(mat1[0])

	mat := NewDmat(dimi, dimk)
	for i := 0; i < dimi; i++ {
		for j := 0; j < dimj; j++ {
			for k := 0; k < dimk; k++ {
				mat[i][k] += mat0[i][j] * mat1[j][k]
			}
		}
	}

	return mat
}

func NewVec(len int) Vec { //Generate a new (zero) vector of length len
	v := make([]float64, len)
	return v
}

func UnitVec(len, dir int) Vec { //Generate a unit vector of length len with dir-th element = 1.0
	v := NewVec(len)
	v[dir] = 1.0
	return v
}

func Ones(len int) Vec { //Generate a vector of ones of length len
	v := NewVec(len)
	for i := range v {
		v[i] = 1.0
	}
	return v
}

func Zeroes(len int) Vec { //Generate a vector of 0 of length len
	return NewVec(len)
}

func multVecVec(vout, v0, v1 Vec) { //element-wise vector multiplication
	for i, v := range v0 {
		vout[i] = v * v1[i]
	}
	return
}

func ScaleVec(vout Vec, s float64, vin Vec) {
	for i, v := range vin {
		vout[i] = s * v
	}
}

func AddVecs(vout, v0, v1 Vec) { //Sum of vectors
	for i := range vout {
		vout[i] = v0[i] + v1[i]
	}
}

func WAddVecs(vout Vec, sca float64, v0, v1 Vec) { //
	for i := range vout {
		vout[i] = sca*v0[i] + v1[i]
	}
}

func DiffVecs(vout, v0, v1 Vec) { //Difference of vectors
	for i := range vout {
		vout[i] = v0[i] - v1[i]
	}
}

func DotVecs(v0, v1 Vec) float64 { //inner product between vectors v0 and v1, use for axis projection
	dot := 0.0
	for i, v := range v0 {
		dot += v * v1[i]
	}
	return dot
}

func Norm2Sq(v Vec) float64 {
	return DotVecs(v, v)
}

func Norm2(v Vec) float64 { //Euclidean Length of vector
	return math.Sqrt(Norm2Sq(v))
}

func NormalizeVec(v Vec) {
	norm := Norm2(v)
	if norm > 0 {
		ScaleVec(v, 1.0/norm, v)
	}
}

func Dist2Vecs(v0, v1 Vec) float64 { //Euclidean distance between 2 vectors squared
	var dev float64
	dist := 0.0
	for i, v := range v0 {
		dev = v - v1[i]
		dist += dev * dev
	}
	return dist
}

func DistVecs(v0, v1 Vec) float64 { //Euclidean distance between 2 vectors
	return math.Sqrt(Dist2Vecs(v0, v1))
}

func DistVecs1(v0, v1 Vec) float64 { //1-norm between 2 vectors
	var dev float64
	dist := 0.0
	for i, v := range v0 {
		dev = v - v1[i]
		dist += math.Abs(dev)
	}
	return dist
}

func Hammingdist(v0, v1 Vec) float64 { //Hamming distance between 2 vectors
	dist := 0.0
	for i, v := range v0 {
		if v != v1[i] {
			dist += 1.0
		}
	}
	return dist
}

func applyFnVec(f func(float64) float64, vec Vec) { //Apply function f to a vector componentwise
	for i, x := range vec {
		vec[i] = f(x)
	}
	return
}

func CopyVec(v Vec) Vec { //makes a copy of a vector
	l := len(v)
	v1 := make(Vec, l)
	copy(v1, v)
	return v1
}

func MinInt(x, y int) int { //Returns minimum of two integers
	if x < y {
		return x
	} else {
		return y
	}
}

func GetMeanVec(vecs []Vec) Vec { // Return the mean vector of array of vectors
	cv := NewVec(len(vecs[0]))
	for _, v := range vecs {
		AddVecs(cv, cv, v)
	}

	fn := 1 / float64(len(vecs))

	ScaleVec(cv, fn, cv)

	return cv
}

func GetVarVec(vecs []Vec) Vec { // Return the variance vector of array of vectors
	lv := len(vecs[0])
	cv := NewVec(lv)
	mv := GetMeanVec(vecs)
	for _, v := range vecs {
		dv := NewVec(lv)
		DiffVecs(dv, v, mv)
		for i, d := range dv {
			cv[i] += d * d
		}
	}

	fn := 1 / float64(len(vecs))

	ScaleVec(cv, fn, cv)

	return cv
}

func SumVec(v Vec) float64 { //Sum of all elements in vector
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

func GetCrossCov(vecs0, vecs1 []Vec, submean0, submean1 bool) (Vec, Vec, Dmat) {
	cv0 := GetMeanVec(vecs0)
	cv1 := GetMeanVec(vecs1)
	ccmat := NewDmat(len(cv0), len(cv1))
	var d0, d1 float64

	for k := range vecs0 {
		for i, c0 := range cv0 {
			if submean0 {
				d0 = vecs0[k][i] - c0
			} else {
				d0 = vecs0[k][i]
			}
			for j, c1 := range cv1 {
				if submean1 {
					d1 = vecs1[k][j] - c1
				} else {
					d1 = vecs1[k][j]
				}
				ccmat[i][j] += d0 * d1
			}
		}
	}

	fn := 1.0 / float64(len(vecs0))
	for i := range cv0 {
		for j := range cv1 {
			ccmat[i][j] *= fn
		}
	}

	return cv0, cv1, ccmat
}

// expect flag = mat.SVDThin or mat.SVDThinU
func GetSVD_opt(ccmat Dmat, flag mat.SVDKind) (*mat.Dense, []float64, *mat.Dense) {
	dim0 := len(ccmat)
	dim1 := len(ccmat[0])
	dim := dim0
	if dim0 > dim1 {
		dim = dim1
	}

	C := mat.NewDense(dim0, dim1, nil)
	for i, ci := range ccmat {
		for j, v := range ci {
			C.Set(i, j, v)
		}
	}

	var svd mat.SVD
	ok := svd.Factorize(C, flag)
	if !ok {
		log.Fatal("SVD failed.")
	}
	U := mat.NewDense(dim0, dim, nil)
	V := mat.NewDense(dim1, dim, nil)
	if flag&mat.SVDThinU == mat.SVDThinU {
		svd.UTo(U)
	}
	if flag&mat.SVDThinV == mat.SVDThinV {
		svd.VTo(V)
	}

	vals := svd.Values(nil)

	return U, vals, V
}

func GetSVD(ccmat Dmat) (*mat.Dense, []float64, *mat.Dense) {
	return GetSVD_opt(ccmat, mat.SVDThin)
}

func ProjectSVD(label string, dirE, dirP *mat.VecDense, data0, data1 Dmat, mean0, mean1 Vec, ccmat Dmat, vals Vec, U, V *mat.Dense) {
	traceP := len(mean0) == len(mean1)
	trace := 0.0
	if traceP {
		for i, v := range ccmat {
			trace += v[i]
		}
	}

	fnorm2 := 0.0
	for _, v := range vals {
		fnorm2 += v * v
	}
	fmt.Printf("%s_FN2,Tr\t%e\t%e\n", label, fnorm2, trace)
	dim0, _ := U.Dims()
	dim1, _ := V.Dims()

	m0 := mat.NewVecDense(dim0, mean0)
	m1 := mat.NewVecDense(dim1, mean1)

	t0 := mat.NewVecDense(dim0, nil)
	t1 := mat.NewVecDense(dim1, nil)
	for i, v := range vals {
		fmt.Printf("%s_vals\t%d\t%e\n", label, i, v)
	}

	fmt.Println("#<YX>=UsV_ali \tcomp\tu.<dY>     \tv.<dX>")
	for i := 0; i < dim0; i++ {
		u := U.ColView(i)
		v := V.ColView(i)
		ue := math.Abs(mat.Dot(dirE, u))
		vp := math.Abs(mat.Dot(dirP, v))
		fmt.Printf("%s_ali\t%d\t%e\t%e\n", label, i, ue, vp)
	}

	fmt.Printf("#<YX>=UsV  \tcomp")
	for i := 0; i < 3; i++ {
		fmt.Printf("\tX.v%-5d\tY.u%-5d", i, i)
	}
	fmt.Printf("\n")
	for k := range data0 {
		fmt.Printf("%s_prj\t%3d", label, k)
		p0 := mat.NewVecDense(dim0, data0[k])
		p1 := mat.NewVecDense(dim1, data1[k])
		for i := 0; i < 3; i++ {
			t0.SubVec(p0, m0)
			u := U.ColView(i)
			y := mat.Dot(t0, u)
			t1.SubVec(p1, m1)
			v := V.ColView(i)
			x := mat.Dot(t1, v)
			if mat.Dot(dirE, u) < 0.0 {
				x *= -1
				y *= -1
			}
			fmt.Printf("\t%e \t%e ", x, y)
		}
		fmt.Printf("\n")
	}
}
