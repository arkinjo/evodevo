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

func PCA2Cue(pop0 *Population, ipca int) Cues {
	pop := pop0.Copy()
	pop.DevPop(0)

	//phenotype in the ancestral environment
	s0 := pop.GetFlatStateVec("P", 0)

	_, _, ccmat := GetCrossCov(s0, s0)
	U, _, _ := GetSVD(ccmat)
	u := U.ColView(ipca)

	dim := u.Len()
	fcues := NewVec(dim)
	for i := 0; i < dim; i++ {
		if math.Signbit(u.AtVec(i)) { //if negative
			fcues[i] = -1
		} else {
			fcues[i] = 1
		}
	}

	cues := NewCues(ncells, nenv)
	clen := ncells + nenv
	for i := 0; i < ncells; i++ {
		copy(cues[i], fcues[i*clen:(i+1)*clen])
	}

	return cues
}

func FlattenEnvs(cues Cues) Vec {
	var v Vec
	for _, cue := range cues {
		v = append(v, cue...)
	}
	return v
}
