package unicell

import (
	"fmt"
	"log"
	"math/rand"
	"os"
)

type Population struct { //Population of individuals
	Env      	Cue //Environment of population in this epoch
	RefEnv		Cue //Environment of population in previous epoch
	Indivs   	[]Indiv
	Fitness  	float64 //Mean fitness of population
	CuePlas		float64 //Mean cue plasticity of population
	ObsPlas		float64 //Mean observed plasticity of population
	Diff		float64 //Mean difference of population
	Utility		float64 //Mean utility of population
}

func NewPopulation(npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	env := NewCue(Nenv) //Initialize environment = zero vector

	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	p := Population{env, env, indivs, 0.0, 0.0, 0.0, 0.0, 0.0}
	return p
}

func (pop *Population) GetMeanFitness() float64 { //average fitness of population
	mf := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mf += indiv.F
	}

	return mf / fn
}

func (pop *Population) GetMeanObsPlasticity() float64 { //average observed plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mp += indiv.Pl
	}

	return mp / fn
}

func (pop *Population) GetMeanCuePlasticity() float64 { //average cue plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mp += indiv.Plc
	}

	return mp / fn
}

func (pop *Population) GetMeanDiff() float64 { //average  of population
	var diff float64
	md := 0.0
	fn := float64(len(pop.Indivs))
	e0 := pop.RefEnv.C //Ancestral environment
	p := NewCue(Nenv)
	for _, indiv := range pop.Indivs {
		p = indiv.Cells[2].P
		diff = dist2Vecs(p.C,e0)
		md += diff
	}

	return md / fn
}


func (pop *Population) GetMeanUtility() float64 { //average utility of population
	mu := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mu += indiv.u
	}

	return mu / fn
}

func (pop *Population) Reproduce(nNewPop int) Population {
	npop := len(pop.Indivs)
	var nindivs []Indiv
	ipop := 0
	cnt := 0
	for ipop < nNewPop && cnt < 1000*nNewPop {
		cnt += 1
		k := rand.Intn(npop)
		l := rand.Intn(npop)
		dad := pop.Indivs[k]
		mom := pop.Indivs[l]
		r0 := rand.Float64()
		r1 := rand.Float64()
		if r0 < dad.F && r1 < mom.F {
			kid0, kid1 := Mate(&dad, &mom)
			//			fmt.Println("kids.G0", len(kid0.G0), len(kid1.G0))

			nindivs = append(nindivs, kid0)
			nindivs = append(nindivs, kid1)
			ipop += 2
		}
	}
	for i := range nindivs {
		nindivs[i].Id = i
	}
	
	new_population := Population{pop.Env, pop.RefEnv, nindivs, 0.0, 0.0, 0.0, 0.0, 0.0} //resets embryonic values to zero!
	/*
		fmt.Println("Reproduce new population: ", ipop)
			for _,indiv := range new_population.Pop {
				fmt.Println("in reproduction", len(indiv.G0))
			}*/

	return new_population
}

func RecEvolve(nstep int, init_pop *Population, epoch int) Population { //Records population fitness and writes file
	var indivenv, refenv Cue
	
	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	//fmt.Fprintln(fout, "Generation", "Population Fitness")
	pop := *init_pop
	ch := make(chan Indiv) //channels for parallelization
	for istep := 1; istep <= nstep; istep++ {
		for _, indiv := range pop.Indivs {
			go func(indiv Indiv) {
				indivenv = pop.Env //Without noise
				refenv = pop.RefEnv
				ch <- indiv.Develop(indivenv, refenv)
			}(indiv)
		}
		for i := range pop.Indivs {
			pop.Indivs[i] = <-ch //Update output results
		}
		pop.Fitness = pop.GetMeanFitness()
		pop.CuePlas = pop.GetMeanCuePlasticity()
		pop.ObsPlas = pop.GetMeanObsPlasticity()
		pop.Utility = pop.GetMeanUtility()
		fmt.Fprintln(fout, epoch, istep, pop.Fitness, pop.CuePlas, pop.ObsPlas, pop.Utility)
		fmt.Println("Evol_step: ", istep, " <Fit>: ", pop.Fitness, "<Pl>:", pop.ObsPlas, "<u>:", pop.Utility) //Prints averages for generation
		pop = pop.Reproduce(MaxPop)
	}
	fout.Close()
	return pop
}
