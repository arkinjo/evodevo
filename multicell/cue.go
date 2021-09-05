package multicell

import (
	//"fmt"
	//"math"
	"math/rand"
)

var DevNoise float64 = 0.00 // Development environment cue noise
var EnvNoise float64 = 0.00 // Selection environment noise


type Cue struct {
	C Vec //Environment cue is a vector
}


type Cues struct {
	Es []Cue //Cue array object
}


func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	tv := NewVec(nenv)         //trait part of vector
	idv := UnitVec(ncells, id) //id part of vector
	v := append(tv, idv...)  //format: cue|id
	cue := Cue{v}
	return cue
}

func (cue *Cue) GetTrait() []float64 { //Extract trait part of cue
	c := cue.C
	tv := c[0:nenv]
	return tv
}

func (cue *Cue) GetId() []float64 { //Extract ID part of cue
	c := cue.C
	idv := c[nenv:] //ID part is appended at end
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
	v := append(tv, idv...)
	cue := Cue{v}

	return cue
}


func (cue *Cue) CopyCue() Cue { //Returns a copy of an environment cue
	c0 := cue.C //this already includes id part
	lc := len(c0)
	cv := make([]float64,lc)
	copy(cv,c0)
	c1 := Cue{cv}

	return c1
}


func (cue *Cue) AddNoise(eta float64) Cue {
	var r float64
	
	tv := cue.GetTrait()
	idv := cue.GetId()

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
	cue1 := Cue{v}

	return cue1
}

func (cue *Cue) ChangeEnv(n int) Cue { // Mutate precisely n bits of environment cue; ignore id part
	env1 := cue.CopyCue() //make a copy of the environmental cue to perform operations without affecting original value
	//Splitting trait part and id part
	tv1 := env1.GetTrait()
	idv := env1.GetId()

	indices := make([]int, len(tv1)) 
	for i := range indices {
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] }) //
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
	
	env1 = Cue{v1}

	return env1
}

func NewCues(ncells, nenv int) Cues {
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = NewCue(nenv, id)
	}
	cs := Cues{vs}
	return cs
}

func RandomEnvs(ncells, nenv int, density float64) Cues { //Randomly generate cue array
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = RandomEnv(nenv, id, density)
	}
	cs := Cues{vs}
	return cs
}

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

func (cues *Cues) AddNoise(eta float64) Cues {
	envs1 := cues.CopyCues()

	for i,c := range envs1.Es {
		envs1.Es[i] = c.AddNoise(eta) //hope this works; I actually don't know what I'm doing here
	}

	return envs1
}

func (cues *Cues) ChangeEnvs(n int) Cues { //Mutates precisely n bits in environment cue vector 'concatenation'
	var ref, cell, cue int
	
	envs1 := cues.CopyCues() //Make a copy to perform operations without changing original value
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
		if envs1.Es[cell].C[cue] == 0 {
			envs1.Es[cell].C[cue] = 1
		} else {
			envs1.Es[cell].C[cue] = 0
		}
	}
	//envs1.Es = cs
	//fmt.Println("Update:",envs1)
	return envs1
}

func (cues *Cues) GetMeanCue() Cue { //elementwise arithmetic mean of environment cue
	//ncells := len(cues)
	//ncues := len(cues)
	cv := NewVec(nenv + ncells)
	for _, env := range cues.Es {
		for j, t := range env.C {
			cv[j] += t / float64(ncells)
		}
	}
	mu := Cue{cv}

	return mu
}

func (cues *Cues) GetCueVar() float64 { //Sum of elementwise variance in environment cue
	mu := cues.GetMeanCue().C
	//ncells := len(cues)
	//ncues := len(cues[0])
	cvec := NewVec(nenv + ncells)
	v := NewVec(nenv + ncells)
	sigma2 := 0.0
	for _, c := range cues.Es {
		cvec = c.C
		diffVecs(v, cvec, mu)
		sigma2 += Veclength2(v) / float64(ncells)
	}
	return sigma2
}
