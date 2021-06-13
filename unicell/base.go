package unicell

import (
	"math"
	"math/rand"
)

var MaxPop int = 200 // population size
var Ngenes int = 100 // number of genes
var Nenv int = 20    // number of environmental cues/phenotypic values per face
var Ncells int = 3 //Three cells, one without cue, one for ancestral environment, one for novel environment

var MaxDevStep int = 200    // Maximum steps for development.
var epsDev float64 = 1.0e-8 // convergence criterion of development.

var GeneLength = 2*Ngenes + 2*Nenv // Length of a gene for Unicellular organism.

var GenomeDensity float64 = 1.0 / float64(Ngenes)

var HalfGenomeDensity float64 = 0.5 * GenomeDensity

var MutationRate float64 = 0.01

var s float64 = 0.25 // selection strength
var Omega float64 = 1.0 // positive parameter of sigmoid, set to limiting to zero (e.g. 1.0e-10) for step function.

var WithCue bool = true // with or without environmental cues. See Develop in indiv.go.

var Traj_Filename string = "test.dat"
var P_Filename string = "phentest.dat" //Expressed phenotypes of population
var G_Filename string = "genetest.dat" //Genome of population

type Spmat = [](map[int]float64) // Sparse matrix is an array of maps.

type Vec = []float64 //Vector is a slice

//var genrand = rand.New(rand.NewSource(99)) //This is bad for concurrency. DO NOT USE!

func SetSeed(seed int64) {
	rand.Seed(seed)
}

func sigmoid(x, omega float64) float64 {
	return 1 / (1 + math.Exp(-x/omega))
}

func sigma(x float64) float64 { //Activation function for development
	return sigmoid(x, Omega)
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

func NewVec(len int) Vec { //Generate a new vector of length len
	v := make([]float64, len)
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

func dist2Vecs(v0, v1 Vec) float64 { //Euclidean distance between 2 vectors squared
	dist := 0.0
	for i, v := range v0 {
		dev := v - v1[i]
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
	nmut := int(MutationRate * float64(nrow*ncol))
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
