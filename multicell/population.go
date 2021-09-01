package multicell

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
//	"io/ioutil" 
)

type Population struct { //Population of individuals
	Gen			int
	Envs      	Cues //Environment of population in this epoch
	RefEnvs		Cues //Environment of population in previous epoch
	Indivs   	[]Indiv
}

func NewPopulation(ncell, npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	envs := NewCues(ncell,Nenv)

	p := Population{0, envs, envs, indivs}
	return p
}

func (pop *Population) GetMeanFitness() float64 { //average fitness of population
	mf := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mf += indiv.Fit
	}

	return mf / fn
}

func (pop *Population) GetMeanObsPlasticity() float64 { //average observed plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mp += indiv.ObsPlas
	}

	return mp / fn
}

func (pop *Population) GetMeanCuePlasticity() float64 { //average cue plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mp += indiv.CuePlas
	}

	return mp / fn
}

func (pop *Population) GetMeanUtility() float64 { //average utility of population
	mu := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mu += indiv.Util
	}

	return mu / fn
}

func (pop *Population) GetMeanPp() float64 { //average degree of polyphenism of population
	mu := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mu += indiv.Pp
	}

	return mu / fn
}

func (pop *Population) GetDiversity() float64 { //To be used after development
	cv := make([]Cue,0)
	for _,ind := range pop.Indivs {
		for _,cell := range ind.Copies[2].Ctypes {
			cv = append(cv, cell.P)
		}
	}
	ccues := Cues{cv}

	return ccues.GetCueVar()
}

func (pop *Population) GetMeanPhenotype(gen int) Cues { //elementwise average phenotype of population; output as slice instead of cue struct
	npop := len(pop.Indivs)
	MeanPhenotype := NewCues(ncells,Nenv)
	pop.DevPop(gen)

	for _,indiv := range(pop.Indivs) {
		for i,c := range(indiv.Copies[2].Ctypes){
			for j,p := range(c.P.C){
				MeanPhenotype.Es[i].C[j] += p/float64(npop)
			}
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
		for i := range Gtilde.Ec {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Ec[i][j] += Gtilde.Ec[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.Eid {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Eid[i][j] += Gtilde.Eid[i][j]/float64(len(pop.Indivs))
			}
		}
		for i := range Gtilde.F {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.F[i][j] += Gtilde.F[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.G {
			for j:=0 ; j<Ngenes ; j++ {
				MeanGenome.G[i][j] += Gtilde.G[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.Hg {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Hg[i][j] += Gtilde.Hg[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.Hh {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Hh[i][j] += Gtilde.Hh[i][j]/float64(len(pop.Indivs))
			} 
		}
		for i := range Gtilde.Pc {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Pc[i][j] += Gtilde.Pc[i][j]/float64(len(pop.Indivs))
			}
		}
		for i := range Gtilde.Pc {
			for j:=0 ; j<Nenv ; j++ {
				MeanGenome.Pid[i][j] += Gtilde.Pid[i][j]/float64(len(pop.Indivs))
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

func (pop *Population) Get_Environment_Axis() Cues { //Choice of axis defined using difference of environment cues
	axlength2 := 0.0

	e := pop.Envs.Es //Cue in novel (present) environment
	e0 := pop.RefEnvs.Es //Cue in ancestral (previous) environment
	v := NewVec(Nenv)
	de := NewCues(ncells,Nenv)

	for i,p := range e {
		diffVecs(v,p.C,e0[i].C)
		axlength2 += Veclength2(v)
		de.Es[i] = Cue{NewVec(ncells),v} //ids must stay the same
	}

	axlength := math.Sqrt(axlength2)
	if axlength == 0 { //if no change in environment cue
		return de
	} else { //normalize
		for i,c := range de.Es{
			for j,p := range c.C {
				de.Es[i].C[j] = p/axlength //normalize to unit vector 
			}
		}
		return de
	}
}

func(pop *Population) Get_Mid_Env() Cues { //Midpoint between ancestral (previous) and novel (current) environment
	e := pop.Envs.Es // novel environment
	e0 := pop.RefEnvs.Es // ancestral environment

	me := NewCues(ncells,Nenv) // midpoint
	
	for i,c := range e {
		for j,v := range c.C {
			me.Es[i].C[j] = (v + e0[i].C[j])/2.0
		}
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
		if r0 < dad.Fit && r1 < mom.Fit {
			kid0, kid1 := Mate(&dad, &mom)
			nindivs = append(nindivs, kid0)
			nindivs = append(nindivs, kid1)
			ipop += 2
		}
	}
	for i := range nindivs {
		nindivs[i].Id = i //Relabels individuals according to position in array
	}
	
	new_population := Population{0, pop.Envs, pop.RefEnvs, nindivs} //resets embryonic values to zero!

	return new_population
}

func (pop *Population) DevPop(gen int) Population {
	novenv := NewCues(ncells,Nenv)
	ancenv := NewCues(ncells,Nenv)

	pop.Gen = gen

	ch := make(chan Indiv) //channels for parallelization
	for _, indiv := range pop.Indivs {
		go func(indiv Indiv) {
			novenv = pop.Envs //Novel environment
			ancenv = pop.RefEnvs //Ancestral environment
			ch <- indiv.CompareDev(novenv, ancenv)
		}(indiv)
	}
	for i := range pop.Indivs {
		pop.Indivs[i] = <-ch //Update output results
	}

	return *pop
}

func Evolve(test bool, tfilename, jsonout, gidfilename string, nstep, epoch int, init_pop *Population) Population { //Records population trajectory and writes files
	var jfilename, id_filename, id, dadid, momid string 
	var Fitness, CuePlas, ObsPlas, Polyp, Div, Util float64
	pop := *init_pop
	//bugfixpop := NewPopulation(len(pop.Envs.Es),len(pop.Indivs))

	if test && gidfilename != ""{ //write genealogy in test mode
		id_filename = fmt.Sprintf("../analysis/%s.dot",gidfilename)
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
				jfilename = fmt.Sprintf("../pops/%s_%d.json",jsonout,pop.Gen)
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
				err = popout.Close()
				if err != nil {
					log.Fatal()
				}
			}

			/*
			if jsonout != "" {//Bugfixing
				jfilename = fmt.Sprintf("../analysis/%s_%d.json",jsonout,pop.Gen)

				popin, err := os.Open(jfilename)
				if err != nil {
					log.Fatal(err)
				}

				byteValue, _ := ioutil.ReadAll(popin)
				err = json.Unmarshal(byteValue, &bugfixpop)
				if err != nil {
					log.Fatal(err)
				}
		
				err = popin.Close()
				if err != nil{
					log.Fatal(err)
				}
				fmt.Println("Nov Env in:",bugfixpop.Envs)
				fmt.Println("Anc Env in:",bugfixpop.RefEnvs)
			}
			*/


		}

		Fitness = pop.GetMeanFitness()
		CuePlas = pop.GetMeanCuePlasticity()
		ObsPlas = pop.GetMeanObsPlasticity()
		Polyp = pop.GetMeanPp()
		Div = pop.GetDiversity()
		Util = pop.GetMeanUtility()

		fout, err := os.OpenFile(tfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}

		fmt.Fprintf(fout,"%d\t%d\t%f\t%e\t%e\t%e\t%e\t%e\n" ,epoch, istep, Fitness, CuePlas, ObsPlas, Polyp, Div, Util)
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf("Evol_step: %d\t <Fit>: %f\t <CPl>:%e\t <OPl>:%e\t <Pp>:%e\t <Div>:%e \t <u>:%e\n ", istep, Fitness, CuePlas, ObsPlas, Polyp, Div, Util)
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
	
	cphen := make(Vec, Nenv)
	mu := pop.Get_Mid_Env()
	Paxis := pop.Get_Environment_Axis() //Bug with something to do with refenv in json file output giving same value to envs and refenvs.
	Projfilename := fmt.Sprintf("../analysis/%s_%d.dat",Filename,gen)

	fout, err := os.OpenFile(Projfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout,"Phenotype\t Genotype")

	for _,indiv := range pop.Indivs {
		pproj, gproj = 0.0, 0.0
		for i,env := range mu.Es { //For each environment cue
			diffVecs(cphen,indiv.Copies[2].Ctypes[i].P.C,env.C) //centralize
			pproj += innerproduct(cphen,Paxis.Es[i].C)
		}
		for i, m := range indiv.Genome.Ec {
			for j, d := range m {
				gproj += d * Gaxis.Ec[i][j]
			}
		}
		for i, m := range indiv.Genome.Eid {
			for j, d := range m {
				gproj += d*Gaxis.Eid[i][j]
			}
		}
		for i, m := range indiv.Genome.F {
			for j, d := range m {
				gproj += d * Gaxis.F[i][j]
			}
		}
		for i, m := range indiv.Genome.G {
			for j, d := range m {
				gproj += d * Gaxis.G[i][j]
			}
		}
		for i, m := range indiv.Genome.Hg {
			for j, d := range m {
				gproj += d * Gaxis.Hg[i][j]
			}
		}
		for i, m := range indiv.Genome.Hh {
			for j, d := range m {
				gproj += d * Gaxis.Hg[i][j]
			}
		}
		for i, m := range indiv.Genome.Pc {
			for j, d := range m {
				gproj += d * Gaxis.Pc[i][j]
			}
		}
		for i, m := range indiv.Genome.Pid {
			for j, d := range m {
				gproj += d * Gaxis.Pid[i][j]
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
