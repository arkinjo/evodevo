package multicell

import (
	//"fmt"
	//"math"
	"math/rand"
)

var DevNoise float64 = 0.00 // Development environment cue noise
var EnvNoise float64 = 0.00 // Selection environment noise

/*
type Cue struct {
	id Vec //id vector of cue
	C Vec //Environment cue is a vector
}
*/
type Cue = Vec

/*
type Cues struct {
	Es []Cue //Cue array object
}
*/

type Cues = []Cue

func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	tv := NewVec(nenv)         //trait part of vector
	idv := UnitVec(ncells, id) //id part of vector
	cue := append(tv, idv...)  //format: cue|id
	return cue
}

func GetTrait(cue Cue) []float64 { //Extract trait part of cue
	tv := cue[0:nenv]
	return tv
}

func GetId(cue Cue) []float64 { //Extract ID part of cue
	idv := cue[nenv:] //ID part is appended at end
	return idv
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
	env := append(tv, idv...)

	return env
}

/* We might not need this function
func CopyCue(c Cue) Cue { //Returns a copy of an environment cue
	c0 := c //this already includes id part
	lc := len(c0)
	cv := make([]float64,lc)
	copy(cv,c0)

	return cv
}
*/
func AddNoisetoCue(cue Cue, eta float64) Cue {
	var r float64
	env1 := make(Cue, nenv)
	copy(env1, cue)
	for i, c := range env1 {
		r = rand.Float64()
		if r < eta {
			if c == 0 {
				env1[i] = 1
			} else {
				env1[i] = 0
			}
		}
	}
	idv := GetId(cue)
	ncue := append(env1, idv...)

	return ncue
}

func ChangeEnv(cue Cue, n int) Cue { // Mutate precisely n bits of environment cue; ignore id part
	//env1 := cue.CopyCue() //make a copy of the environmental cue to perform operations without affecting original value
	//v1 := env1.C

	v1 := make([]float64, nenv)
	copy(v1, cue)

	indices := make([]int, len(v1))
	for i := range indices {
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] }) //
	mutindices := make([]int, n)
	for i := range mutindices {
		mutindices[i] = indices[i]
	}
	for _, i := range mutindices {
		if v1[i] == 0 {
			v1[i] = 1
		} else {
			v1[i] = 0
		}
	}

	idv := GetId(cue)
	env1 := append(v1, idv...)

	return env1
}

func NewCues(ncells, nenv int) Cues {
	cs := make([]Cue, ncells)
	for id := range cs {
		cs[id] = NewCue(nenv, id)
	}
	//cs := Cues{vs}
	return cs
}

func RandomEnvs(ncells, nenv int, density float64) Cues { //Randomly generate cue array
	cs := make([]Cue, ncells)
	for id := range cs {
		cs[id] = RandomEnv(nenv, id, density)
	}
	//cs := Cues{vs}
	return cs
}

/*
func (cues *Cues) CopyCues() Cues{
	ncells := len(cues.Es)
	vs := make([]Cue,ncells)
	//copy(vs,cues.Es) //Original
	for i,c := range cues.Es { //'Proactive copying'
		vs[i] = c.CopyCue() //No bugs with this implementation here
	}
	cs := Cues{vs}
	return cs
}
*/

func AddNoisetoCues(cues Cues, eta float64) Cues {
	//envs1 := cues.CopyCues()

	envs1 := make([]Cue, ncells)
	for i := range envs1 {
		envs1[i] = AddNoisetoCue(cues[i], eta)
	}

	return envs1
}

func ChangeEnvs(cues Cues, n int) Cues { //Mutates precisely n bits in environment cue vector 'concatenation'
	var ref, cell, cue int
	//envs1 := cues.CopyCues() //Make a copy to perform operations without changing original value
	envs1 := make([]Cue, ncells)

	copy(envs1, cues)
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
	v := NewVec(nenv + ncells)
	sigma2 := 0.0
	for _, cvec := range cues {
		diffVecs(v, cvec, mu)
		sigma2 += Veclength2(v) / float64(ncells)
	}
	return sigma2
}
