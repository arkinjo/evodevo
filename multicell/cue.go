package multicell

import (
	"math/rand"
)

var rand_cue = rand.New(rand.NewSource(99))

type Cue = Vec //Environment cue is a special kind of vector

type Cues = []Cue //Cue array object

func SetSeedCue(seed int64) {
	rand_cue.Seed(seed)
}

func RandomEnv(density float64) Cue { //Fake up a boolean environment vector for celltype id
	v := NewVec(nenv)
	for i := range v {
		if rand_cue.Float64() < density {
			v[i] = cueMag
		} else {
			v[i] = -cueMag
		}
	}
	return v
}

// Mutate precisely n bits of environment cue; ignore id part
func ChangeEnv(cue Cue, n int) Cue {
	env1 := CopyVec(cue)
	if n == 0 {
		return env1
	}
	indices := make([]int, nenv)
	for i := range indices {
		indices[i] = i
	}
	rand_cue.Shuffle(len(indices), func(i, j int) { indices[i], indices[j] = indices[j], indices[i] })
	for _, i := range indices[0:n] {
		env1[i] = -env1[i]
	}

	return env1
}

func NewCues(ncells, nenv int) Cues {
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = NewVec(nenv)
	}
	return vs
}

//Randomly generate cue array
func RandomEnvs(ncells int, density float64) Cues {
	vs := make([]Cue, ncells)
	for id := range vs {
		vs[id] = RandomEnv(density)
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

func AddNoise2CueNormal(cue_out, cue Cue, eta float64) {
	copy(cue_out, cue)
	for i, t := range cue {
		// Don't use rand_cue here. Use the system rand instead.
		cue_out[i] = t + eta*rand.NormFloat64()
	}

	return
}

func AddNoise2CueFlip(cue_out, cue Cue, eta float64) {
	copy(cue_out, cue)
	for i, t := range cue {
		// Don't use rand_cue here. Use the system rand instead.
		if rand.Float64() < eta {
			cue_out[i] = -t
		}
	}

	return
}

func ChangeEnvs(cues Cues, n int) Cues { //Flips precisely n bits in each environment cue
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

func GetSelEnvs(cues Cues) Cues {
	cues1 := make([]Cue, len(cues))
	for i, cue := range cues {
		cues1[i] = NewVec(nsel)
		copy(cues1[i], cue[0:nsel])
	}

	return cues1
}

func FlattenEnvs(cues Cues) Vec {
	var v Vec
	for _, cue := range cues {
		v = append(v, cue...)
	}
	return v
}
