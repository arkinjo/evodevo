package multicell

import (
	//"fmt"
	//"math"
	"math/rand"
)

var DevNoise float64 = 0.00 // Development environment cue noise
var EnvNoise float64 = 0.00 // Selection environment noise

type Cue = Vec //Environment cue is a vector

type Cues = []Cue //Cue array object

func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	tv := NewVec(nenv)         //trait part of vector
	idv := UnitVec(ncells, id) //id part of vector
	v := append(tv, idv...)    //format: cue|id
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

func GetId(cue Cue) int { //Extract id number
	var id int
	idv := GetIdVec(cue)
	for i, v := range idv {
		if v != 0 {
			id = i
		}
	}
	return id

}

func RandomEnv(nenv, id int, density float64) Cue { //Fake up a boolean environment vector for celltype id
	var r float64

	tv := make([]float64, nenv)
	for i := range tv {
		r = rand.Float64()
		if r < density {
			tv[i] = 1
		} else {
			tv[i] = 0
		}
	}
	idv := UnitVec(ncells, id)
	v := append(tv, idv...)

	return v
}

/* This is now redundant after introducing CopyVec function
func CopyCue(cue Cue) Cue { //Returns a copy of an environment cue
	c0 := cue //this already includes id part
	lc := len(c0)
	cv := make([]float64, lc)
	copy(cv, c0)

	return cv
}
*/

func AddNoisetoCue(cue Cue, eta float64) Cue {
	var r float64

	tv := GetTrait(cue)
	idv := GetIdVec(cue)

	for i, t := range tv {
		r = rand.Float64()
		if r < eta {
			if t == 0 {
				tv[i] = 1
			} else {
				tv[i] = 0
			}
		}
	}
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
		if tv1[i] == 0 {
			tv1[i] = 1
		} else {
			tv1[i] = 0
		}
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

func AddNoisetoCues(cues Cues, eta float64) Cues {
	envs1 := CopyCues(cues)

	for i, c := range envs1 {
		envs1[i] = AddNoisetoCue(c, eta) //hope this works; I actually don't know what I'm doing here
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
		if envs1[cell][cue] == 0 {
			envs1[cell][cue] = 1
		} else {
			envs1[cell][cue] = 0
		}
	}
	//envs1.Es = cs
	//fmt.Println("Update:",envs1)
	return envs1
}

func GetMeanCue(cues Cues) Cue { //elementwise arithmetic mean of environment cue
	//ncells := len(cues)
	//ncues := len(cues)
	cv := NewVec(nenv + ncells)
	for _, env := range cues {
		for j, t := range env {
			cv[j] += t / float64(ncells)
		}
	}

	return cv
}

func GetCueVar(cues Cues) float64 { //Sum of elementwise variance in environment cue
	mu := GetMeanCue(cues)
	//ncells := len(cues)
	//ncues := len(cues[0])
	//cvec := NewVec(nenv + ncells)
	v := NewVec(nenv + ncells)
	sigma2 := 0.0
	for _, c := range cues {
		//cvec = c //possibly redundant line
		diffVecs(v, c, mu)
		sigma2 += Veclength2(v) / float64(ncells)
	}
	return sigma2
}
