package multicell

import (
	"math"
	"math/rand"
)

type Settings struct {
	MaxPop  int     // Maximum number of individuals in population
	NCells  int     // Number of cell types
	ELayer  bool    // e present?
	FLayer  bool    // f present?
	HLayer  bool    // h present?
	JLayer  bool    //  J present?
	Pfback  bool    // P feedback to E layer
	SDNoise float64 // stdev of environmental noise
}

var default_settings = Settings{1000, 1, true, true, true, true, true, 0.05}

var maxPop int = 1000 // population size
var ngenes int = 200  // number of genes
var nenv int = 40     // number of environmental cues/phenotypic values per face
var ncells int = 1    //number of cell types/phenotypes to be trained simultaneously; not exported
var pheno_feedback bool = false
var devNoise float64 = 0.05

const maxDevStep int = 200    // Maximum steps for development.
const ccStep int = 5          // Number of steady steps for convergence
const epsDev float64 = 1.0e-5 // Convergence criterion of development.
const eps float64 = 1.0e-50

const alphaEMA = 2.0 / (1.0 + 5.0) // exponential moving average/variance

var fullGeneLength = 4*ngenes + 2*nenv + 2*ncells // Length of a gene for Unicellular organism.
var genelength int                                //calculated from layers present or absent.

const funcspergene float64 = 1.0 //average number of functions per gene
var GenomeDensity float64 = funcspergene / float64(ngenes)
var CueResponseDensity float64 = -math.Log(eps) / float64(ngenes)

var HalfGenomeDensity float64 = 0.5 * GenomeDensity

const baseMutationRate float64 = 0.01 // default probability of mutation of genome
var mutRate float64                   //declaration
const baseSelStrength float64 = 20    // default selection strength; to be normalized by number of cells
const selDevStep float64 = 20.0       // Developmental steps for selection

//var selStrength float64             //declaration; Selection strength per unit cue
//var minFitness float64
const minWagnerFitness float64 = 0.01

var Omega float64 = 1.0 // positive parameter of sigmoid, set to limiting to zero (e.g. 1.0e-10) for step function.

var withCue bool = false // with or without environmental cues.
var cuestrength float64  // cue strength
var epig bool = true     // Epigenetic marker layer
var hoc bool = true      // Higher order complexes layer
var hoi bool = false
var jstrength float64 //declaration

// theoretical standard deviation of matrix elements
var sdE float64
var sdF float64
var sdG float64
var sdH float64
var sdJ float64
var sdP float64

//Remark: defaults to full model!

type Vec = []float64 //Vector is a slice
type Dmat = []Vec

//var genrand = rand.New(rand.NewSource(99)) //This is bad for concurrency. DO NOT USE!

func SetSeed(seed int64) {
	rand.Seed(seed)
}

func SetParams(s Settings) { //Define whether each layer or interaction is present in model
	maxPop = s.MaxPop
	withCue = s.ELayer
	if s.ELayer {
		cuestrength = 1.0
	} else {
		cuestrength = 0.0
	}

	if s.JLayer {
		jstrength = 1.0
	} else {
		jstrength = 0.0
	}

	epig = s.FLayer
	hoc = s.HLayer
	hoi = s.JLayer
	pheno_feedback = s.Pfback
	ncells = s.NCells
	devNoise = s.SDNoise

	genelength = ngenes + (nenv + ncells) //G and P layers present by default
	if s.ELayer {
		genelength += nenv + ncells
	}

	if s.FLayer {
		genelength += ngenes
	}
	if s.HLayer {
		genelength += ngenes
		if s.JLayer {
			genelength += ngenes
		}
	}

	mutRate = baseMutationRate * float64(fullGeneLength) / float64(genelength) //to compensate for layer removal.

	//initializing theoretical standard deviations for entries of each matrix
	sdE = math.Sqrt(cuestrength / (CueResponseDensity * float64(nenv+ncells) * (1 + cuestrength)))
	if s.Pfback {
		sdE *= math.Sqrt(0.5) // rescale to account for feedback
	}

	sdG = 1 / math.Sqrt(GenomeDensity*float64(ngenes)*(1+cuestrength))

	sdF = math.Sqrt(math.Pi / (float64(ngenes) * GenomeDensity))
	sdH = 1 / math.Sqrt(GenomeDensity*float64(ngenes)*(1+jstrength))
	sdJ = math.Sqrt(jstrength / (GenomeDensity * float64(ngenes) * (1 + jstrength)))
	sdP = 1 / math.Sqrt(CueResponseDensity*float64(ngenes))
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

func scale(x float64) float64 {
	return cuestrength * x
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
	return 6.0 * math.Atan(x/1.73205080756887729352744634150587236694) / math.Pi
}

/*
func scaledatan(x, omega float64) float64 {
	return 2.0 * math.Atan(omega*x) / math.Pi
}
*/

func relu(x, omega float64) float64 {
	if x < 0 {
		return 0
	} else {
		return omega * x
	}
}

func sigmaf(x float64) float64 { //Activation function for epigenetic markers
	return lecunatan(x)
}

func sigmag(x float64) float64 { //Activation function for gene expression levels
	return lecunatan(x)
}

func sigmah(x float64) float64 { //Activation function for higher order complexes
	return lecunatan(x) //abstract level of amount of higher order complexes
}

func rho(x float64) float64 { //Function for converting gene expression into phenotype
	return lecunatan(x)
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

func DiffVecs(vout, v0, v1 Vec) { //Difference of vectors
	for i := range vout {
		vout[i] = v0[i] - v1[i]
	}
}

func Innerproduct(v0, v1 Vec) float64 { //inner product between vectors v0 and v1, use for axis projection
	dot := 0.0
	for i, v := range v0 {
		dot += v * v1[i]
	}
	return dot
}

func Norm2Sq(v Vec) float64 {
	return Innerproduct(v, v)
}

func Norm2(v Vec) float64 { //Euclidean Length of vector
	return math.Sqrt(Norm2Sq(v))
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

func GetCrossCov(vecs0, vecs1 []Vec) (Vec, Vec, Dmat) {
	cv0 := GetMeanVec(vecs0)
	cv1 := GetMeanVec(vecs1)
	ccmat := NewDmat(len(cv0), len(cv1))

	for k := range vecs0 {
		for i, c0 := range cv0 {
			d0 := vecs0[k][i] - c0
			for j, c1 := range cv1 {
				d1 := vecs1[k][j] - c1
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
