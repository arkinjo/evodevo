package multicell


import (
	//"fmt"
	//"math"
	"math/rand"
)

var DevNoise float64 = 0.00 // Development environment cue noise
var EnvNoise float64 = 0.00 // Selection environment noise


type Cue struct {
	id Vec //id vector of cue
	C Vec //Environment cue is a vector
}

type Cues struct {
	Es []Cue //Cue array object
}

func NewCue(nenv, id int) Cue { //Initialize a new cue type object
	cv := NewVec(nenv)
	idv := UnitVec(ncells, id)
	cue := Cue{idv,cv}
	return cue
}


func RandomEnv(nenv, id int, density float64) Cue { //Fake up a boolean environment vector
	var r float64
	env := NewCue(nenv, id)
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
	lc := len(c0)
	cv := make([]float64,lc)
	copy(cv,c0)

	id0 := cue.id
	lid := len(id0)
	idv := make([]float64,lid)
	copy(idv,id0)

	c1 := Cue{idv,cv}
	return c1
}

func (cue *Cue) AddNoise(eta float64) Cue {
	var r float64
	env1 := cue.CopyCue()
	v1 := env1.C //ignore id vector!
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


func (cue *Cue) ChangeEnv(n int) Cue {// Mutate precisely n bits of environment cue; ignore id part
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
	for id := range vs { 
		vs[id] = NewCue(id,nenv)
	}
	cs := Cues{vs}
	return cs
}

func RandomEnvs(ncells, nenv int, density float64) Cues { //Randomly generate cue array
	vs := make([]Cue,ncells)
	for id := range vs { 
		vs[id] = RandomEnv(nenv,id,density)
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
	for i := range envs1.Es{
		envs1.Es[i] = envs1.Es[i].AddNoise(eta)
	}
	return envs1
}

func (cues *Cues) ChangeEnvs(n int) Cues { //Mutates precisely n bits in environment cue vector 'concatenation'
	var ref, cell, cue int
	//fmt.Println("Input:",cues)
	envs1 := cues.CopyCues() //Make a copy to perform operations without changing original value
	//fmt.Println("Copy:",envs1)
	cs := envs1.Es
	N := Nenv*ncells
	indices := make([]int,N)
	for i:= range indices {
		indices[i] = i
	}

	rand.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] })
	mutcells := make([]int,0)
	mutcues := make([]int,0)
	for i := 0; i < n; i++ {
		ref = indices[i]
		cell = ref/Nenv //integer division
		cue = ref%Nenv
		mutcells = append(mutcells,cell)
		mutcues = append(mutcues,cue)
	}
	//fmt.Println("Cell index:",mutcells)
	//fmt.Println("Cue index:",mutcues)
	for j,cell := range mutcells {
		cue = mutcues[j]
		if cs[cell].C[cue] == 0{
			cs[cell].C[cue] = 1
		} else {
			cs[cell].C[cue] = 0
		}
	}
	envs1.Es = cs
	//fmt.Println("Update:",envs1)
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