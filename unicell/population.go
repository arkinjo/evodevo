package unicell

import (
	"fmt"
	"log"
	"math/rand"
	"os"
)

type Population struct { //Population of individuals
	Gen			int
	Env      	Cue //Environment of population in this epoch
	RefEnv		Cue //Environment of population in previous epoch
	Indivs   	[]Indiv
}

func NewPopulation(npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	env := NewCue(Nenv) //Initialize environment = zero vector

	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	p := Population{0, env, env, indivs}
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

func (pop *Population) Reproduce(nNewPop int) Population { //Makes new generation of individuals
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
		nindivs[i].Id = i //Relabels individuals according to position in array
	}
	
	new_population := Population{0, pop.Env, pop.RefEnv, nindivs} //resets embryonic values to zero!

	return new_population
}

func (pop *Population) DevPop(gen int) Population {
	var indivenv, refenv Cue

	pop.Gen = gen

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

func Evolve(test bool, tfilename, pfilename, gfilename string, nstep, epoch int, init_pop *Population) Population { //Records population fitness and writes file
	var filename string
	var Gtilde Genome 
	var Fitness, CuePlas, ObsPlas, Util float64
	pop := *init_pop

	for istep := 1; istep <= nstep; istep++ {		
		pop.DevPop(istep)
		if test { //Dump phenotypes and genotypes in test mode
			if pfilename != "" {
				filename = fmt.Sprintf("%s%d_%d.dat",pfilename,epoch,istep)
				fout, err := os.OpenFile(filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
				if err != nil {
					log.Fatal(err)
				}
				for _, indiv := range pop.Indivs {
					for _,trait := range indiv.Cells[2].P.C {
						fmt.Fprintf(fout, "%e\t", trait )
					}
					fmt.Fprint(fout,"\n")
		
				}
				err = fout.Close()
				if err != nil {
					log.Fatal(err)
				}
			}
			if gfilename != "" {
				filename = fmt.Sprintf("%s%d_%d.dat",gfilename,epoch,istep)
				fout, err := os.OpenFile(filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
				if err != nil {
					log.Fatal(err)
				}
				for _, indiv := range pop.Indivs {
					Gtilde = indiv.Genome
					for i := range Gtilde.G {
						for j:=0 ; j < Ngenes ; j++ {
							fmt.Fprintf(fout, "%e\t",Gtilde.G[i][j])
						}
					}
					Gtilde = indiv.Genome
					for i := range Gtilde.E {
						for j:=0 ; j < Nenv ; j++ {
							fmt.Fprintf(fout, "%e\t",Gtilde.E[i][j])
						}
					}
					Gtilde = indiv.Genome
					for i := range Gtilde.P {
						for j:=0 ; j < Nenv ; j++ {
							fmt.Fprintf(fout, "%e\t",Gtilde.P[i][j])
						}
					}
					Gtilde = indiv.Genome
					for i := range Gtilde.Z {
						for j:=0 ; j < Ngenes ; j++ {
							fmt.Fprintf(fout, "%e\t",Gtilde.Z[i][j])
						}
					}
					fmt.Fprint(fout,"\n")
				}
				err = fout.Close()
				if err != nil {
					log.Fatal(err)
				}
			}
		}

		Fitness = pop.GetMeanFitness()
		CuePlas = pop.GetMeanCuePlasticity()
		ObsPlas = pop.GetMeanObsPlasticity()
		Util = pop.GetMeanUtility()

		fout, err := os.OpenFile(tfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}

		fmt.Fprintf(fout,"%d\t%d\t%f\t%e\t%e\t%e\n" ,epoch, istep, Fitness, CuePlas, ObsPlas, Util)
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf("Evol_step: %d\t <Fit>: %f\t <Epg>:%e\t <Pl>:%e\t <u>:%e\n ", istep, Fitness, CuePlas, ObsPlas, Util)
		pop = pop.Reproduce(MaxPop)
	}
	return pop
}

func (pop *Population) DOT_Ids(DotFilename string, gen int) { //Dump in DOT language friendly format
	var id, dadid, momid string
	Filename := fmt.Sprintf("%s.dot",DotFilename)
	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout,"") //Insert first line here, Digraph etc.
	for _, indiv := range pop.Indivs {
		id = fmt.Sprintf("g%d:id%d",gen,indiv.Id)
		dadid = fmt.Sprintf("g%d:id%d",gen-1,indiv.DadId) //Dad and mom from previous generation
		momid = fmt.Sprintf("g%d:id%d",gen-1,indiv.MomId)
		fmt.Fprintf(fout,"%s -> {%s %s}\n",id,dadid,momid) //Use child -> parent convention
		}
	fmt.Fprintln(fout,"") //Insert final line here, end of graph generating function
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

func (pop *Population) Dump_Phenotypes(Filename string, gen int) { //Extracts phenotypes from population
	
	pop.DevPop(gen)

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	for _, indiv := range pop.Indivs {
		for _,trait := range indiv.Cells[2].P.C {
			fmt.Fprintf(fout, "%e\t", trait )
		}
		fmt.Fprint(fout,"\n")
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

func (pop *Population) Dump_Genotypes(Filename string) { //Extracts genomes from population
	var Gtilde Genome

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}

	for _, indiv := range pop.Indivs {
		Gtilde = indiv.Genome
		for i := range Gtilde.G {
			for j:=0 ; j < Ngenes ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.G[i][j])
			}
		}
		Gtilde = indiv.Genome
		for i := range Gtilde.E {
			for j:=0 ; j < Nenv ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.E[i][j])
			}
		}
		Gtilde = indiv.Genome
		for i := range Gtilde.P {
			for j:=0 ; j < Nenv ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.P[i][j])
			}
		}
		Gtilde = indiv.Genome
		for i := range Gtilde.Z {
			for j:=0 ; j < Ngenes ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.Z[i][j])
			}
		}
		fmt.Fprint(fout,"\n")
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

