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
	Plasticity	float64 //Mean plasticity of population
	Utility		float64 //Mean utility of population
}

func NewPopulation(npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	//env := RandomEnv(Nenv, 0.5)
	env := NewCue(Nenv) //Initialize environment = zero vector

	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	p := Population{env, env, indivs, 0.0, 0.0, 0.0}
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

func (pop *Population) GetIndivFitness() []float64 { //fitness of each individual in population
	var Flist []float64
	//Flist := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		Flist = append(Flist, indiv.F)
	}
	return Flist
}

func (pop *Population) GetMeanPlasticity() float64 { //average plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mp += indiv.Pl
	}

	return mp / fn
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
	
	new_population := Population{pop.Env, pop.RefEnv, nindivs, 0.0, 0.0, 0.0} //resets embryonic values to zero!
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
		pop.Plasticity = pop.GetMeanPlasticity()
		pop.Utility = pop.GetMeanUtility()
		fmt.Fprintln(fout, epoch, istep, pop.Fitness, pop.Plasticity, pop.Utility)
		fmt.Println("Evol_step: ", istep, " <Fit>: ", pop.Fitness, "<Pl>:", pop.Plasticity, "<u>:", pop.Utility) //Prints averages for generation
		pop = pop.Reproduce(MaxPop)
	}
	fout.Close()
	return pop
}

/*
func RecIndivEvolve(nstep int, init_pop *Population, epoch int) Population { //Records population fitness and individual fitness and writes file
	// Warning! This function DOES NOT work!
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
				//indivenv := pop.Env //Without noise
				indivenv := AddNoise(pop.Env, NoisyEta) //Environment of individual may be different
				ch <- indiv.Develop(indivenv)
			}(indiv)
		}
		for i := range pop.Indivs {
			pop.Indivs[i] = <-ch //Update output results
		}
		pop.Fitness = pop.GetMeanFitness()
		//Fslice := pop.GetIndivFitness()
		//fmt.Println(istep,":",len(Fslice),"=",MaxPop,"is",len(Fslice)==MaxPop) //Diagnostics; Parallel processing seems to create issues with this strategy.
		fmt.Fprint(fout, epoch, istep, pop.Fitness)
		for _,indiv := range pop.Indivs{
			fmt.Fprint(fout,indiv.Fitness)
		}
		fmt.Fprint(fout,"\r\n")
		fmt.Println("Evol_step: ", istep, " <fitness>: ", pop.Fitness) //Prints average fitness for generation
		pop = pop.Reproduce(MaxPop)
	}
	err = fout.Close()
	if err != nil{
		log.Fatal(err)
	}
	return pop
}
*/
