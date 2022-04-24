package multicell

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
)

var rand_cue = rand.New(rand.NewSource(99))

var devNoise float64 = 0.05 // Development environment cue noise
//var envNoise float64 = 0.00 // Selection environment noise

type Cue = Vec //Environment cue is a special kind of vector

type Cues = []Cue //Cue array object

var ZeroEnvs Cues

func init() {
	ZeroEnvs = make([]Cue, ncells)
	for i := range ZeroEnvs {
		tv := make([]float64, nenv)
		idv := UnitVec(ncells, i)
		ZeroEnvs[i] = append(tv, idv...)
	}
}

func SetSeedCue(seed int64) {
	rand_cue.Seed(seed)
}

func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	tv := NewVec(nenv)         //trait part of vector
	idv := UnitVec(ncells, id) //id part of vector
	for j, u := range idv {
		if u != 1 {
			idv[j] = -1 //pm 1 representation for id as well
		}
	}
	v := append(tv, idv...) //format: cue|id
	return v
}

func GetTrait(cue Cue) []float64 { //Extract trait part of cue
	tv := cue[0:nenv]

	return tv
}

func GetIdVec(cue Cue) []float64 { //Extract ID part of cue
	idv := cue[nenv:] //ID part is appended at end
	return idv
}

func GetId(cue Cue) int { //Extract cue id number
	var id int
	idv := GetIdVec(cue)
	for i, v := range idv {
		if v == 1 {
			id = i
			break //exit loop prematurely after hitting id
		}
	}
	return id

}

func RandomEnv(nenv, id int, density float64) Cue { //Fake up a boolean environment vector for celltype id
	var r float64

	tv := make([]float64, nenv)
	for i := range tv {
		r = rand_cue.Float64()
		if r < density { //density is probability of 1
			tv[i] = 1
		} else {
			tv[i] = -1
		}
	}
	idv := UnitVec(ncells, id)
	for j, u := range idv {
		if u != 1 {
			idv[j] = -1 //pm 1 representation for id as well
		}
	}
	v := append(tv, idv...)

	return v
}

func SetNoise(eta float64) {
	devNoise = eta
}

func ChangeEnv(cue Cue, n int) Cue { // Mutate precisely n bits of environment cue; ignore id part
	env1 := CopyVec(cue) //make a copy of the environmental cue to perform operations without affecting original value
	//Splitting trait part and id part
	tv1 := GetTrait(env1)
	idv := GetIdVec(env1)

	indices := make([]int, len(tv1))
	for i := range indices {
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] })
	mutindices := make([]int, n)
	for i := range mutindices {
		mutindices[i] = indices[i]
	}
	for _, i := range mutindices {
		tv1[i] = -tv1[i] //Flip bit
	}

	v1 := append(tv1, idv...) //glue back id part

	return v1
}

func NewCues(ncells, nenv int) Cues {
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = NewCue(nenv, id)
	}
	return vs
}

func RandomEnvs(ncells, nenv int, density float64) Cues { //Randomly generate cue array
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = RandomEnv(nenv, id, density)
	}
	return vs
}

func CopyCues(cues Cues) Cues {
	ncells := len(cues)
	vs := make([]Cue, ncells)
	for i, c := range cues { //'Proactive copying'
		vs[i] = CopyVec(c)
	}
	return vs
}

func AddNoise2Cue(cue Cue, eta float64) Cue {
	tv := GetTrait(cue)
	idv := GetIdVec(cue)
	v := CopyVec(tv)
	for i, t := range tv {
		v[i] = t + eta*rand.NormFloat64()
	}
	v = append(v, idv...)

	return v
}

func AddNoise2Cues(cues Cues, eta float64) Cues {
	envs1 := CopyCues(cues)
	for i, c := range envs1 {
		envs1[i] = AddNoise2Cue(c, eta) //hope this works; I actually don't know what I'm doing here
	}

	return envs1
}

func ChangeEnvs(cues Cues, n int) Cues { //Mutates precisely n bits in environment cue vector 'concatenation'
	var ref, cell, cue int

	envs1 := CopyCues(cues) //Make a copy to perform operations without changing original value
	N := nenv * ncells
	indices := make([]int, N)
	for i := range indices {
		indices[i] = i
	}

	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] })
	mutcells := make([]int, 0)
	mutcues := make([]int, 0)
	for i := 0; i < n; i++ {
		ref = indices[i]
		cell = ref / nenv //integer division; evaluates to cell index
		cue = ref % nenv  //modular arithmetic; evaluates to trait index; is always < nenv
		mutcells = append(mutcells, cell)
		mutcues = append(mutcues, cue)
	}
	//fmt.Println("Cell index:",mutcells)
	//fmt.Println("Cue index:",mutcues)
	for j, cell := range mutcells {
		cue = mutcues[j]
		if envs1[cell][cue] == -1 {
			envs1[cell][cue] = 1
		} else {
			envs1[cell][cue] = -1
		}
	}
	//envs1.Es = cs
	//fmt.Println("Update:",envs1)
	return envs1
}

func ChangeEnvs2(cues Cues, n int) Cues { //Flips precisely n bits in each environment cue
	cues1 := CopyCues(cues)
	for i, cue := range cues {
		cues1[i] = ChangeEnv(cue, n)
	}
	return cues1
}

func GetMeanCue(cues Cues) Cue { //elementwise arithmetic mean of environment cue
	cv := NewVec(nenv + ncells)
	for _, env := range cues {
		for j, t := range env {
			cv[j] += t
		}
	}

	fn := 1 / float64(len(cues))

	for j, t := range cv {
		cv[j] = t * fn
	}
	return cv
}

func GetCueVar(cues Cues) float64 { //Sum of elementwise variance in environment cue
	mu := GetMeanCue(cues)
	v := NewVec(nenv + ncells)
	sigma2 := 0.0
	for _, c := range cues {
		diffVecs(v, c, mu)
		sigma2 += Veclength2(v)
	}

	sigma2 /= float64(len(cues))

	return sigma2
}

func PCAtoCue(pcafilename string) ([]Cues, [][][]float64) { //Converts PCA vectors containing only trait part into environment cue vectors
	var prdir, cell, trait, r int
	var x, y float64
	var str string

	PCAenvs := make([]Cues, nenv*ncells) //Number of PCAs = Length of trait concatenation
	for i := range PCAenvs {
		PCAenvs[i] = NewCues(ncells, nenv)
	}

	PCAvecs := make([][][]float64, nenv*ncells*ncells)
	for i := range PCAvecs {
		PCAvecs[i] = make([][]float64, nenv*ncells)
		for j := range PCAvecs[i] {
			PCAvecs[i][j] = make([]float64, nenv)
		}
	}

	traits := make([]float64, 0) //nenv and ncells are fixed and given; this is the concatenation of trait parts of cue
	values := make([]float64, 0)

	filename := fmt.Sprintf("../analysis/%s.dat", pcafilename)
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	scanner := bufio.NewScanner(file)

	for scanner.Scan() { //error in scanner?!
		str = scanner.Text()
		if x, err = strconv.ParseFloat(str, 64); err == nil {
			values = append(values, x)
			if math.Signbit(x) { //if negative
				y = -1.0
			} else {
				y = 1.0
			}
			traits = append(traits, y)
		}
	}
	fmt.Println(len(traits) == len(values))
	/*
		if len(traits) != nenv*ncells { //Is this even needed?
			err := errors.New("input incompatable with size")
			return PCAenvs, err
		}
	*/
	for i, t := range traits {
		prdir = i / (nenv * ncells) //take advantage of integer division
		r = i % (nenv * ncells)     //remainder
		cell = r / nenv
		trait = r % nenv
		//fmt.Printf("index:%d\tprdir:%d\tr:%d\tcindex:%d\ttindex:%d\ttrait:%f\n", i, prdir, r, cell, trait, t)
		PCAenvs[prdir][cell][trait] = t
		PCAvecs[prdir][cell][trait] = values[i]
	}

	return PCAenvs, PCAvecs //, nil
}
