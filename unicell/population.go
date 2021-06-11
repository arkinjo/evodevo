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
	Utility		float64 //Mean utility of population
}

func NewPopulation(npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	env := NewCue(Nenv) //Initialize environment = zero vector

	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	p := Population{env, env, indivs, 0.0, 0.0, 0.0, 0.0}
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
			nindivs = append(nindivs, kid0)
			nindivs = append(nindivs, kid1)
			ipop += 2
		}
	}
	for i := range nindivs {
		nindivs[i].Id = i
	}
	
	new_population := Population{pop.Env, pop.RefEnv, nindivs, 0.0, 0.0, 0.0, 0.0} //resets embryonic values to zero!

	return new_population
}

func (pop *Population) DevPop() Population {
	var indivenv, refenv Cue

	ch := make(chan Indiv) //channels for parallelization
	for _, indiv := range pop.Indivs {
		go func(indiv Indiv) {
			indivenv = pop.Env //Novel environment
			refenv = pop.RefEnv //Ancestral environment
			ch <- indiv.Develop(indivenv, refenv)
		}(indiv)
	}
	for i := range pop.Indivs {
		pop.Indivs[i] = <-ch //Update output results
	}

	return *pop
}



func RecEvolve(nstep, epoch int, init_pop *Population) Population { //Records population fitness and writes file
	
	fout, err := os.OpenFile(Traj_Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}

	pop := *init_pop

	for istep := 1; istep <= nstep; istep++ {
		
		pop.DevPop()

		pop.Fitness = pop.GetMeanFitness()
		pop.CuePlas = pop.GetMeanCuePlasticity()
		pop.ObsPlas = pop.GetMeanObsPlasticity()
		pop.Utility = pop.GetMeanUtility()
		fmt.Fprintln(fout, epoch, istep, pop.Fitness, pop.CuePlas, pop.ObsPlas, pop.Utility)
		fmt.Println("Evol_step: ", istep, " <Fit>: ", pop.Fitness, "<Epg>:", pop.CuePlas , "<Pl>:", pop.ObsPlas, "<u>:", pop.Utility) //Prints averages for generation
		pop = pop.Reproduce(MaxPop)
	}
	fout.Close()
	return pop
}

func (pop *Population) Get_Phenotypes(Filename string) { //Extracts phenotypes from population
	var str_trait string
	
	pop.DevPop()

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	for _, indiv := range pop.Indivs {
		//fmt.Fprintln(fout, indiv.Cells[2].P.C)
		for _,trait := range indiv.Cells[2].P.C {
			str_trait = fmt.Sprint(trait)
			fmt.Fprintf(fout, str_trait + "\t" )
		}
		fmt.Fprint(fout,"\n")

	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}


func (pop *Population) Get_Genotypes(Filename string) { //Extracts genomes from population
	var str_val string
	var Gtilde Genome

	pop.DevPop()

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	for _, indiv := range pop.Indivs {
		Gtilde = indiv.Genome 
		for i, row := range Gtilde.G {
			for j := range row {
				str_val = fmt.Sprint(Gtilde.G[i][j])
				fmt.Fprintf(fout, str_val + "\t" )
			}
		}
		for i, row := range Gtilde.E {
			for j := range row {
				str_val = fmt.Sprint(Gtilde.E[i][j])
				fmt.Fprintf(fout, str_val + "\t" )
			}
		}
		for i, row := range Gtilde.P {
			for j := range row {
				str_val = fmt.Sprint(Gtilde.P[i][j])
				fmt.Fprintf(fout, str_val + "\t" )
			}
		}
		for i, row := range Gtilde.Z {
			for j := range row {
				str_val = fmt.Sprint(Gtilde.Z[i][j])
				fmt.Fprintf(fout, str_val + "\t" )
			}
		}
		fmt.Fprint(fout,"\n")
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

