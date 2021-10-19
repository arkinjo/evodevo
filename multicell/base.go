package multicell

import (
	"math"
	"math/rand"
)

var maxPop int = 1000 // population size
var ngenes int = 200  // number of genes
var nenv int = 5      // number of environmental cues/phenotypic values per face
var ncells int = 1    //number of cell types/phenotypes to be trained simultaneously; not exported

const maxDevStep int = 1000   // Maximum steps for development.
const ccStep int = 5          //number of steady steps for convergence
const epsDev float64 = 1.0e-6 // convergence criterion of development.

var fullGeneLength = 4*ngenes + 2*nenv + 2*ncells // Length of a gene for Unicellular organism.
var genelength int                                //calculated from layers present or absent.

const funcspergene float64 = 3.0 //average number of functions per gene
var GenomeDensity float64 = funcspergene / float64(ngenes)
var CueResponseDensity float64 = -math.Log(epsDev) / float64(ngenes)

var HalfGenomeDensity float64 = 0.5 * GenomeDensity

const baseMutationRate float64 = 0.01 // default probability of mutation of genome
var mutRate float64                   //declaration
const baseSelStrength float64 = 0.25  // default selection strength; to be normalized by number of cells
var selStrength float64               //declaration; Remark: not used in L1 norm fitness.
var Omega float64 = 1.0               // positive parameter of sigmoid, set to limiting to zero (e.g. 1.0e-10) for step function.

var withCue bool = false // with or without environmental cues.
var cuestrength float64  //declaration
var epig bool = true     // Epigenetic marker layer
var hoc bool = true      // Higher order complexes layer
var hoi bool = false
var hoistrength float64 //declaration

var sde float64  //declaration
var sdf float64  //declaration
var sdg float64  //declaration
var sdhg float64 //declaration
var sdhh float64 //declaration
var sdp float64  //declaration

//Remark: defaults to full model!

type Spmat struct {
	Ncol int                 // number of columns
	Mat  [](map[int]float64) // Sparse matrix is an array of maps.
}

type Vec = []float64 //Vector is a slice

//var genrand = rand.New(rand.NewSource(99)) //This is bad for concurrency. DO NOT USE!

func SetSeed(seed int64) {
	rand.Seed(seed)
}

func SetMaxPop(n int) {
	maxPop = n
}

func SetNcells(n int) {
	ncells = n
	selStrength = baseSelStrength / float64(n)
}

func SetLayers(ce, ch float64, epigm, HOC bool) { //Define whether each layer or interaction is present in model
	cuestrength = ce //strength of environment cue
	hoistrength = ch //strength of interactions between higher order complexes

	//withCue = cue //Whether environment cue has effect on development
	epig = epigm //Layer representing epigenetic markers
	hoc = HOC    //Layer representing higher-order complexes
	//hoi = HOI    //Allow interaction between higher-order complexes

	genelength = ngenes + nenv + ncells //G and P layers present by default
	if ce != 0 {
		withCue = true
		genelength += nenv + ncells
	} else {
		withCue = false
	}

	if epig {
		genelength += ngenes
	}
	if hoc {
		genelength += ngenes
		if ch != 0 {
			hoi = true
			genelength += ngenes
		} else {
			hoi = false
		}
	}

	mutRate = baseMutationRate * float64(fullGeneLength) / float64(genelength) //to compensate for layer removal.

	//initializing theoretical standard deviations
	sdg = 1 / math.Sqrt(GenomeDensity*float64(ngenes)*(1+cuestrength))
	sde = math.Sqrt(cuestrength / (CueResponseDensity * float64(nenv+ncells) * (1 + cuestrength)))
	sdf = math.Sqrt(math.Pi / (float64(ngenes) * GenomeDensity))
	sdhg = 1 / math.Sqrt(GenomeDensity*float64(ngenes)*(1+hoistrength))
	sdhh = math.Sqrt(hoistrength / (GenomeDensity * float64(ngenes) * (1 + hoistrength)))
	sdp = 1 / math.Sqrt(CueResponseDensity*float64(ngenes))
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
	return relu(x, 1.0)
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
				sp.Mat[i][j] = rand.NormFloat64() / sd //Normalize wrt sd of theoretical output
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

func NewVec(len int) Vec { //Generate a new (zero) vector of length len
	v := make([]float64, len)
	return v
}

func UnitVec(len, dir int) Vec { //Generate a unit vector of length len with dir-th element = 1.0
	v := NewVec(len)
	v[dir] = 1.0
	return v
}

func multMatVec(vout Vec, mat Spmat, vin Vec) { //Matrix multiplication
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

func multMatVec_T(vout Vec, mat Spmat, vin Vec) { //Matrix transposition and then multiplication
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

func addVecs(vout, v0, v1 Vec) { //Sum of vectors
	for i := range vout {
		vout[i] = v0[i] + v1[i]
	}
}

func diffVecs(vout, v0, v1 Vec) { //Difference of vectors
	for i := range vout {
		vout[i] = v0[i] - v1[i]
	}
}

func innerproduct(v0, v1 Vec) float64 { //inner product between vectors v0 and v1, use for axis projection
	dot := 0.0
	for i, v := range v0 {
		dot += v * v1[i]
	}
	return dot
}

func Veclength2(v Vec) float64 {
	return innerproduct(v, v)
}

func Veclength(v Vec) float64 { //Euclidean Length of vector
	return math.Sqrt(Veclength2(v))
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

func (mat *Spmat) mutateSpmat(density, sd float64) { //mutating a sparse matrix
	nrow := len(mat.Mat)
	nmut := int(mutRate * float64(nrow*mat.Ncol))
	for n := 0; n < nmut; n++ {
		i := rand.Intn(nrow)
		j := rand.Intn(mat.Ncol)
		r := rand.Float64()
		delete(mat.Mat[i], j)
		if r < density {
			mat.Mat[i][j] = rand.NormFloat64() / sd //normalize wrt theoretical sd of output
		}
	}
	//Note: This implementation has non-zero probability of choosing same element to be mutated twice.
	return
}

func CopyVec(v Vec) Vec { //makes a copy of a vector
	l := len(v)
	v1 := make(Vec, l)
	copy(v1, v)
	return v1
}

