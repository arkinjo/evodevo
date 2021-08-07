package unicell


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

func NewCue(nenv int) Cue { //Initialize a new cue type object
	cv := NewVec(nenv)
	cue := Cue{cv}
	return cue
}


func RandomEnv(nenv int, density float64) Cue { //Fake up a boolean environment vector
	var r float64
	env := NewCue(nenv)
	v := make([]float64, nenv)
	for i := range v {
		r = rand.Float64()
		if r < density {
			v[i] = 1
		} else {
			v[i] = 0
		}
	}
	env.C = v

	return env
}

func (cue *Cue) CopyCue() Cue { //Returns a copy of an environment cue
	c0 := cue.C
	l := len(c0)
	v := make([]float64,l)
	copy(v,c0)
	c1 := Cue{v}
	return c1
}

func (cue *Cue) AddNoise(eta float64) Cue {
	var r float64
	env1 := cue.CopyCue()
	v1 := env1.C
	for i, c := range v1 {
		r = rand.Float64()
		if r < eta {
			if c == 0 {
				v1[i] = 1
			} else {
				v1[i] = 0
			}
		}
	}
	env1.C = v1

	return env1
}


func (cue *Cue) ChangeEnv(n int) Cue {// Mutate precisely n bits of environment cue.
	env1 := cue.CopyCue() //make a copy of the environmental cue to perform operations without affecting original value
	v1 := env1.C
	indices := make([]int,len(v1))
	for i := range indices{
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] }) //
	mutindices := make([]int,n)
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

	env1.C = v1
	
	return env1
}

func NewCues(ncells, nenv int) Cues {
	vs := make([]Cue,ncells)
	for i := range vs { 
		vs[i] = NewCue(nenv)
	}
	cs := Cues{vs}
	return cs
}

func RandomEnvs(ncells, nenv int, density float64) Cues { //Randomly generate cue array
	vs := make([]Cue,ncells)
	for i := range vs { 
		vs[i] = RandomEnv(nenv,density)
	}
	cs := Cues{vs}
	return cs
}

func (cues *Cues) CopyCues() Cues{
	ncells := len(cues.Es)
	vs := make([]Cue,ncells)
	copy(vs,cues.Es)
	cs := Cues{vs}
	return cs
}

func (cues *Cues) AddNoise(eta float64) Cues {
	envs1 := cues.CopyCues()
	for i := range envs1.Es{
		envs1.Es[i] = envs1.Es[i].AddNoise(eta)
	}
	return envs1
}

func (cues *Cues) ChangeEnv(n int) Cues { //Mutates precisely n bits in environment cue vector 'concatenation'
	var ref, cell, cue int
	envs1 := cues.CopyCues()
	cs := envs1.Es
	N := Nenv*Ncells
	indices := make([]int,N)
	for i:= range indices{
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] }) //
	mutcells := make([]int,0)
	mutcues := make([]int,0)
	for i := 0; i<n; i++ {
		ref = indices[i]
		cell = ref/Nenv
		cue = ref%Nenv
		mutcells = append(mutcells,cell)
		mutcues = append(mutcues,cue)
	}
	for j,cell := range mutcells {
		cue = mutcues[j]
		if cs[cell].C[cue] == 0{
			cs[cell].C[cue] = 1
		} else {
			cs[cell].C[cue] = 0
		}
	}
	return envs1
}


func (cues *Cues) GetMeanCue() Vec { //Mean of environment cue
	ncells := len(cues.Es)
	ncues := len(cues.Es[0].C)
	cv := NewVec(ncues)
	for _,env := range cues.Es {
		for j,c := range env.C{
			cv[j] += c/float64(ncells)
		}
	}
	return cv
}

func (cues *Cues) GetCueVar() float64 { //Variance in environment cue
	mu := cues.GetMeanCue()
	ncells := len(cues.Es)
	ncues := len(cues.Es[0].C)
	v := NewVec(ncues)
	sigma2 := 0.0
	for _,cvec := range cues.Es {
		diffVecs(v,cvec.C,mu)
		sigma2 += Veclength2(v)/float64(ncells)
	}
	return sigma2
}