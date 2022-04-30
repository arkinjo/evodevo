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

type Cue = Vec //Environment cue is a special kind of vector

type Cues = []Cue //Cue array object

func SetSeedCue(seed int64) {
	rand_cue.Seed(seed)
}

func NewIDVec(n, id int) Vec {
	v := NewVec(n)
	for i := range v {
		if i == id {
			v[i] = 1
		} else {
			v[i] = -1
		}
	}
	return v
}
func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	tv := NewVec(nenv)          //trait part of vector
	idv := NewIDVec(ncells, id) //id part of vector
	v := append(tv, idv...)     //format: cue|id
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
	idv := NewIDVec(ncells, id)
	v := append(tv, idv...)

	return v
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

func AddNoise2Cue(cue_out, cue Cue, eta float64) {
	tv := GetTrait(cue)
	for i, t := range tv {
		cue_out[i] = t + eta*rand.NormFloat64()
	}

	return
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

func GetCueVar(cues Cues) float64 { //Sum of elementwise variance in environment cue
	mu := GetMeanVec(cues)
	v := NewVec(nenv + ncells)
	sigma2 := 0.0
	for _, c := range cues {
		DiffVecs(v, c, mu)
		sigma2 += Norm2Sq(v)
	}

	sigma2 /= float64(len(cues))

	return sigma2
}

//Converts PCA vectors containing only trait part into environment cue vectors
func PCAtoCue(pcafilename string) ([]Cues, [][][]float64) {
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

func FlattenEnvs(cues Cues) Vec {
	var v Vec
	for _, cue := range cues {
		v = append(v, cue...)
	}
	return v
}
