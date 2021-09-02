package multicell

import (
	"math"
	"math/rand"
)

var MaxPop int = 1000 // population size
var Ngenes int = 400  // number of genes
var Nenv int = 20     // number of environmental cues/phenotypic values per face
var ncells int = 1    //number of cell types/phenotypes to be trained simultaneously; not exported

var MaxDevStep int = 200    // Maximum steps for development.
var epsDev float64 = 1.0e-8 // convergence criterion of development.

var fullGeneLength = 5*Ngenes + 2*Nenv + 2*ncells // Length of a gene for Unicellular organism.

var GenomeDensity float64 = 1.0 / float64(Ngenes)

var HalfGenomeDensity float64 = 0.5 * GenomeDensity

const baseMutationRate float64 = 0.01 // default probability of mutation of genome
var mutRate float64                   //declaration
const baseSelStrength float64 = 0.25  // default selection strength; to be normalized by number of cells
var selStrength float64               //declaration
var Omega float64 = 1.0               // positive parameter of sigmoid, set to limiting to zero (e.g. 1.0e-10) for step function.

var withCue bool = true // with or without environmental cues.
var epig bool = true    // Epigenetic marker layer
var hoc bool = true     // Higher order complexes layer
var hoi bool = true     // Interaction between higher order complexes

//Remark: defaults to full model!

type Spmat = [](map[int]float64) // Sparse matrix is an array of maps.

type Vec = []float64 //Vector is a slice

//var genrand = rand.New(rand.NewSource(99)) //This is bad for concurrency. DO NOT USE!

func SetSeed(seed int64) {
	rand.Seed(seed)
}

func SetNcells(n int) {
	ncells = n
	selStrength = baseSelStrength / float64(n)
}

func SetLayers(cue, epigm, HOC, HOI bool) { //Define whether each layer or interaction is present in model
	withCue = cue //Whether environment cue has effect on development
	epig = epigm  //Layer representing epigenetic markers
	hoc = HOC     //Layer representing higher-order complexes
	hoi = HOI     //Allow interaction between higher-order complexes

	genelength := 2*Ngenes + Nenv + ncells
	if cue {
		genelength += Nenv + ncells
	}
	if epig {
		genelength += Ngenes
	}
	if hoc {
		genelength += Ngenes
		if hoi {
			genelength += Ngenes
		}
	}

	mutRate = baseMutationRate * float64(fullGeneLength) / float64(genelength) //to compensate for layer removal.
}

func GetNcells() int {
	return ncells
}

func sigmoid(x, omega float64) float64 {
	return 1 / (1 + math.Exp(-x/omega))
}

func relu(x, omega float64) float64 {
	if x < 0 {
		return 0
	} else {
		return omega * x
	}
}

func sigmaf(x float64) float64 { //Activation function for epigenetic markers
	return sigmoid(x, Omega)
}

func sigmag(x float64) float64 { //Activation function for gene expression levels
	return relu(x, Omega)
}

func sigmah(x float64) float64 { //Activation function for higher order complexes
	if hoi { // prevent explosion by bounding
		return sigmoid(x, Omega)
	} else {
		return relu(x, Omega)
	}
}

func rho(x float64) float64 { //Function for converting gene expression into phenotype
	return sigmoid(x, Omega)
}

func NewSpmat(nrow, ncol int, density float64) Spmat { //Generate a new sparse matrix
	mat := make([](map[int]float64), nrow)
	for i := range mat {
		mat[i] = make(map[int]float64)
		for j := 0; j < ncol; j++ {
			r := rand.Float64()
			if r < density {
				mat[i][j] = rand.NormFloat64()
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

func multMatVec(vout Vec, mat Spmat, vin Vec) { //Matrix multiplication
	for i := range vout {
		vout[i] = 0.0
	}

	for i, m := range mat {
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
	for i, m := range mat {
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

func dist2Vecs(v0, v1 Vec) float64 { //Euclidean distance between 2 vectors squared
	var dev float64
	dist := 0.0
	for i, v := range v0 {
		dev = v - v1[i]
		dist += dev * dev
	}
	return dist
}

func distVecs(v0, v1 Vec) float64 { //Euclidean distance between 2 vectors
	return math.Sqrt(dist2Vecs(v0, v1))
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

func mutateSpmat(mat Spmat, ncol int) { //mutating a sparse matrix
	nrow := len(mat)
	nmut := int(baseMutationRate * float64(nrow*ncol))
	for n := 0; n < nmut; n++ {
		i := rand.Intn(nrow)
		j := rand.Intn(ncol)
		r := rand.Float64()
		delete(mat[i], j)
		if r < GenomeDensity {
			mat[i][j] = rand.NormFloat64()
		}
	}
	return
}
