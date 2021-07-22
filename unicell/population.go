package unicell

import (
	"fmt"
	"log"
	"math/rand"
	"encoding/json"
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

func (pop *Population) GetMeanPhenotype(gen int) Vec { //elementwise average phenotype of population
	npop := len(pop.Indivs)
	MeanPhenotype := make(Vec, npop)
	pop.DevPop(gen)

	for _,indiv := range(pop.Indivs) {
		for i,p := range(indiv.Cells[2].P.C){
			MeanPhenotype[i] += p/float64(npop)
		}
	}
	return MeanPhenotype
}

//func (pop *Population) CentralizePhenotypes(gen int, Centre Vec)

func (pop *Population) GetMeanGenome() Genome { //elementwise average genome of population
	var Gtilde Genome
	
	MeanGenome := NewGenome()
	for _,indiv := range pop.Indivs {
		Gtilde = indiv.Genome
		for i := range Gtilde.G {
			for j:=0 ; j<Ngenes ; j++ {
				MeanGenome.G[i][j] += Gtilde.G[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.E {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.E[i][j] += Gtilde.E[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.P {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.P[i][j] += Gtilde.P[i][j]/float64(len(pop.Indivs))
			}
		}
		for i := range Gtilde.Z {
			for j:=0 ; j<Ngenes ; j++ {
				MeanGenome.Z[i][j] += Gtilde.Z[i][j]/float64(len(pop.Indivs))
			} 
		}
	}
	return MeanGenome
}

func (pop *Population) Get_Environment_Axis() Vec { //Choice of axis defined using difference of environment cues
	e := pop.Env.C
	e0 := pop.RefEnv.C
	de := NewVec(len(e))
	diffVecs(de, e, e0)
	difflength := Veclength(de)
	for i := range de { //normalize
		de[i] = de[i]/difflength 
	}
	return de
}

func(pop *Population) Get_Mid_Env() Vec { //Midpoint between environments
	e := pop.Env.C
	e0 := pop.RefEnv.C
	me := NewVec(len(e))
	addVecs(me, e, e0)
	for i := range me {
		me[i] = me[i]/2
	}
	return me
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

func Evolve(test bool, tfilename, jsonout, gidfilename string, nstep, epoch int, init_pop *Population) Population { //Records population fitness and writes file
	var id_filename, id, dadid, momid string 
	var Fitness, CuePlas, ObsPlas, Util float64
	pop := *init_pop

	if test && gidfilename != ""{
		id_filename = fmt.Sprintf("%s.dot",gidfilename)
		fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(fout,"digraph G {")
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}
	}

	for istep := 1; istep <= nstep; istep++ {		
		pop.DevPop(istep)
		if test{
			if gidfilename!= "" && istep>1 { //Genealogy not defined for first generation
				fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
				if err != nil {
					log.Fatal(err)
				}
				for _, indiv := range pop.Indivs {
					id = fmt.Sprintf("g%d:id%d",pop.Gen,indiv.Id)
					dadid = fmt.Sprintf("g%d:id%d",pop.Gen-1,indiv.DadId) //Dad and mom from previous generation
					momid = fmt.Sprintf("g%d:id%d",pop.Gen-1,indiv.MomId)
					fmt.Fprintf(fout,"\t%s -> {%s, %s}\n",id,dadid,momid) //Use child -> parent convention
				}
			}
			if jsonout!="" { //Export JSON population of each generation in test mode
				jfilename := fmt.Sprintf("%s_%d.json",jsonout,pop.Gen)
				jsonpop, err := json.Marshal(pop) //JSON encoding of population as byte array
				if err != nil {
					log.Fatal(err)
				}
				popout, err := os.OpenFile(jfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create json file
				if err != nil {
					log.Fatal(err)
				}
				_, err = popout.Write(jsonpop)
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
	if test && gidfilename!= ""{
		fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(fout,"}")
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}
	}
	return pop
}

func (pop *Population) Dump_Projections(Filename string, gen int, Gaxis Genome) {
	var pproj, gproj float64
	pop.DevPop(gen)

	cphens := make([]Vec,len(pop.Indivs))
	mu := pop.Get_Mid_Env()
	for i,indiv := range(pop.Indivs){
		diffVecs(cphens[i], indiv.Cells[2].P.C, mu)//centralize
	}
	Paxis := pop.Get_Environment_Axis()

	Projfilename := fmt.Sprintf("%s_%d.dat",Filename,gen)

	fout, err := os.OpenFile(Projfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout,"Phenotype\t Genotype")

	for i,indiv := range(pop.Indivs){
		pproj = innerproduct(cphens[i],Paxis)
		gproj = 0.0 //initialize
		for i, m := range indiv.Genome.G {
			for j, d := range m {
				gproj += d * Gaxis.G[i][j]
			}
		}
		for i, m := range indiv.Genome.E {
			for j, d := range m {
				gproj += d * Gaxis.E[i][j]
			}
		}
		for i, m := range indiv.Genome.P {
			for j, d := range m {
				gproj += d * Gaxis.P[i][j]
			}
		}
		for i, m := range indiv.Genome.Z {
			for j, d := range m {
				gproj += d * Gaxis.Z[i][j]
			}
		}
		fmt.Fprintf(fout,"%e\t %e\n",pproj,gproj)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

/* Code beyond this point may be considered as obsolete.
func (pop *Population) Project_Phenotypes(Filename string, gen int, axis Vec) []float64 {
	var proj float64
	pop.DevPop(gen)
	
	projs := make([]float64,len(pop.Indivs))

	for _,indiv := range(pop.Indivs){
		proj = innerproduct(indiv.Cells[2].P.C,axis)
		projs = append(projs, proj)
		
	}

	return projs
}

func (pop *Population) Project_Genomes(Filename string, axis Genome) []float64 {
	var proj float64
	projs := make([]float64,len(pop.Indivs))

	for _,indiv := range(pop.Indivs){
		proj = 0.0 //initialize
		for i, m := range indiv.Genome.G {
			for j, d := range m {
				proj += d * axis.G[i][j]
			}
		}
		for i, m := range indiv.Genome.E {
			for j, d := range m {
				proj += d * axis.E[i][j]
			}
		}
		for i, m := range indiv.Genome.P {
			for j, d := range m {
				proj += d * axis.P[i][j]
			}
		}
		for i, m := range indiv.Genome.Z {
			for j, d := range m {
				proj += d * axis.Z[i][j]
			}
		}
		projs = append(projs, proj)
	}
	return projs
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
		for i := range Gtilde.E {
			for j:=0 ; j < Nenv ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.E[i][j])
			}
		}
		for i := range Gtilde.P {
			for j:=0 ; j < Nenv ; j++ {
				fmt.Fprintf(fout, "%e\t",Gtilde.P[i][j])
			}
		}
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
*/
