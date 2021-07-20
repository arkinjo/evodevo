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

func NewCue(Length int) Cue { //Initialize a new cue type object
	cv := NewVec(Length)
	cue := Cue{cv}
	return cue
}


func RandomEnv(n int, density float64) Cue { //Fake up a boolean environment vector
	var r float64
	env := NewCue(Nenv)
	v := make([]float64, n)
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

