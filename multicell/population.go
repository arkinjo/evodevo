package multicell

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	//	"io/ioutil"
)

type Population struct { //Population of individuals
	Gen     int
	Envs    Cues //Environment of population in this epoch
	RefEnvs Cues //Environment of population in previous epoch
	Indivs  []Indiv
}

func NewPopulation(ncell, npop int) Population { //Initialize new population
	indivs := make([]Indiv, npop)
	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

	envs0 := NewCues(ncell, nenv)
	envs1 := NewCues(ncell, nenv)

	p := Population{0, envs0, envs1, indivs}
	return p
}

func (pop *Population) SetWagnerFitness() { //compute normalized fitness value similar to Wagner (1996).
	var mf float64
	for _, indiv := range pop.Indivs {
		if mf < indiv.Fit {
			mf = indiv.Fit
		}
	}
	for _, indiv := range pop.Indivs {
		indiv.WagFit = indiv.Fit / mf
	}
}

func (pop *Population) RandomizeGenome() {
	for _, indiv := range pop.Indivs { //Sets genome of every individual to zero
		indiv.Genome.Randomize()
	}
}

func (pop *Population) ClearGenome() {
	for _, indiv := range pop.Indivs { //Sets genome of every individual to zero
		indiv.Genome.Clear()
	}
}

func (pop *Population) Copy() Population {
	npop := len(pop.Indivs)
	ncell := len(pop.Indivs[0].Copy().Copies[0].Ctypes) //number of cells
	//fmt.Println("Copying ",ncell,"-cell individuals")
	pop1 := NewPopulation(ncell, npop)
	pop1.Gen = pop.Gen
	pop1.Envs = CopyCues(pop.Envs)
	pop1.RefEnvs = CopyCues(pop.RefEnvs)
	for i, indiv := range pop.Indivs {
		pop1.Indivs[i] = indiv.Copy()
	}
	return pop1
}

func (pop *Population) SortPopIndivs() {
	sort.Slice(pop.Indivs, func(i, j int) bool { return pop.Indivs[i].Id < pop.Indivs[j].Id }) //Hopefully this works
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
	cv := make([]Cue, 0)
	for _, ind := range pop.Indivs {
		for _, cell := range ind.Copies[2].Ctypes {
			cv = append(cv, cell.P)
		}
	}
	div := GetCueVar(cv)

	return div
}

func (pop *Population) GetMeanPhenotype(gen int) Cues { //elementwise average phenotype of population; output as slice instead of cue struct
	npop := len(pop.Indivs)
	MeanPhenotype := NewCues(ncells, nenv)
	pop.DevPop(gen)

	for _, indiv := range pop.Indivs {
		for i, c := range indiv.Copies[2].Ctypes {
			for j, p := range c.P {
				MeanPhenotype[i][j] += p / float64(npop)
			}
		}
	}
	return MeanPhenotype
}

//func (pop *Population) CentralizePhenotypes(gen int, Centre Vec)

func (pop *Population) GetMeanGenome() Genome { //elementwise average genome of population
	var Gtilde Genome
	fnpop := 1.0 / float64(len(pop.Indivs))
	MeanGenome := NewGenome()

	for _, indiv := range pop.Indivs {
		Gtilde = indiv.Genome
		if withCue {
			for i, m := range Gtilde.E.Mat {
				for j, v := range m {
					MeanGenome.E.Mat[i][j] += v * fnpop
				}
			}
		}
		if epig {
			for i, m := range Gtilde.F.Mat {
				for j, v := range m {
					MeanGenome.F.Mat[i][j] += v * fnpop
				}
			}
		}
		for i, m := range Gtilde.G.Mat {
			for j, v := range m {
				MeanGenome.G.Mat[i][j] += v * fnpop
			}
		}
		if hoc {
			for i, m := range Gtilde.Hg.Mat {
				for j, v := range m {
					MeanGenome.Hg.Mat[i][j] += v * fnpop
				}
			}
			if hoi {
				for i, m := range Gtilde.Hh.Mat {
					for j, v := range m {
						MeanGenome.Hh.Mat[i][j] += v * fnpop
					}
				}
			}
		}
		for i, m := range Gtilde.P.Mat {
			for j, v := range m {
				MeanGenome.P.Mat[i][j] += v * fnpop
			}
		}
		/*
			for i, m := range Gtilde.Z.Mat {
				for j, v := range m {
					MeanGenome.Z.Mat[i][j] += v * fnpop
				}
			}
		*/
	}

	return MeanGenome
}

func (pop *Population) Get_Environment_Axis() Cues { //Choice of axis defined using difference of environment cues
	axlength2 := 0.0

	e := pop.Envs     //Cue in novel (present) environment
	e0 := pop.RefEnvs //Cue in ancestral (previous) environment
	v := NewVec(nenv + ncells)
	de := NewCues(ncells, nenv)

	for i, p := range e {
		diffVecs(v, p, e0[i])
		axlength2 += Veclength2(v)
		de[i] = v //ids must stay the same
	}

	axlength := math.Sqrt(axlength2)
	if axlength == 0 { //if no change in environment cue
		return de
	} else { //normalize
		for i, c := range de {
			for j, p := range c {
				de[i][j] = p / axlength //normalize to unit vector
			}
		}
		return de
	}
}

func (pop *Population) Get_Mid_Env() Cues { //Midpoint between ancestral (previous) and novel (current) environment
	e := pop.Envs     // novel environment
	e0 := pop.RefEnvs // ancestral environment

	me := NewCues(ncells, nenv) // midpoint

	for i, c := range e {
		for j, v := range c {
			me[i][j] = (v + e0[i][j]) / 2.0
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
		if r0 < dad.WagFit && r1 < mom.WagFit {
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
	novenv := NewCues(ncells, nenv)
	ancenv := NewCues(ncells, nenv)

	pop.Gen = gen

	ch := make(chan Indiv) //channels for parallelization
	for _, indiv := range pop.Indivs {
		go func(indiv Indiv) {
			novenv = pop.Envs    //Novel environment
			ancenv = pop.RefEnvs //Ancestral environment
			ch <- indiv.CompareDev(novenv, ancenv)
		}(indiv)
	}
	for i := range pop.Indivs {
		pop.Indivs[i] = <-ch //Update output results
	}
	//We might need a sorter here.

	return *pop
}

func Evolve(test bool, tfilename, jsonout, gidfilename string, nstep, epoch int, init_pop *Population) Population { //Records population trajectory and writes files
	var jfilename, id_filename, id, dadid, momid string
	var Fitness, CuePlas, ObsPlas, Polyp, Div, Util float64
	pop := *init_pop
	//bugfixpop := NewPopulation(len(pop.Envs.Es),len(pop.Indivs))

	if test && gidfilename != "" { //write genealogy in test mode
		id_filename = fmt.Sprintf("../analysis/%s.dot", gidfilename)
		fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(fout, "digraph G {")
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}
	}

	for istep := 1; istep <= nstep; istep++ {
		pop.DevPop(istep)
		if test {
			if gidfilename != "" && istep > 1 { //Genealogy not defined for first generation
				fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
				if err != nil {
					log.Fatal(err)
				}
				for _, indiv := range pop.Indivs {
					id = fmt.Sprintf("g%d:id%d", pop.Gen, indiv.Id)
					dadid = fmt.Sprintf("g%d:id%d", pop.Gen-1, indiv.DadId) //Dad and mom from previous generation
					momid = fmt.Sprintf("g%d:id%d", pop.Gen-1, indiv.MomId)
					fmt.Fprintf(fout, "\t%s -> {%s, %s}\n", id, dadid, momid) //Use child -> parent convention
				}
			}
			if jsonout != "" { //Export JSON population of each generation in test mode
				jfilename = fmt.Sprintf("../pops/%s_%d.json", jsonout, pop.Gen)
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

		fmt.Fprintf(fout, "%d\t%d\t%f\t%e\t%e\t%e\t%e\t%e\n", epoch, istep, Fitness, CuePlas, ObsPlas, Polyp, Div, Util)
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf("Evol_step: %d\t <Fit>: %f\t <CPl>:%e\t <OPl>:%e\t <Pp>:%e\t <Div>:%e \t <u>:%e\n ", istep, Fitness, CuePlas, ObsPlas, Polyp, Div, Util)
		pop = pop.Reproduce(maxPop)
	}
	if test && gidfilename != "" {
		fout, err := os.OpenFile(id_filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(fout, "}")
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}
	}
	return pop
}

func (pop *Population) Dump_Projections(Filename string, gen int, Gaxis Genome, Paxis Cues) {
	var ancpproj, novpproj, gproj float64
	pop.DevPop(gen) //Not needed for bugfixing

	anccphen := make(Vec, nenv+ncells)
	novcphen := make(Vec, nenv+ncells)
	mu := pop.Get_Mid_Env()
	//fmt.Println("Middle environment : ", mu)
	Projfilename := fmt.Sprintf("../analysis/%s_%d.dat", Filename, gen)

	fout, err := os.OpenFile(Projfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "NovPhen \t AncPhen \t Genotype")

	for _, indiv := range pop.Indivs {
		ancpproj, novpproj, gproj = 0.0, 0.0, 0.0

		for i, env := range mu { //For each environment cue
			//copy(indiv.Copies[2].Ctypes[i].P, pop.Envs[i])       //Bugfixing test
			diffVecs(novcphen, indiv.Copies[2].Ctypes[i].P, env) //centralize
			novpproj += innerproduct(novcphen, Paxis[i])
		}
		for i, env := range mu { //For each environment cue
			//copy(indiv.Copies[1].Ctypes[i].P, pop.RefEnvs[i])    //Bugfixing test
			diffVecs(anccphen, indiv.Copies[1].Ctypes[i].P, env) //centralize
			ancpproj += innerproduct(anccphen, Paxis[i])         //Plot phenotype when pulled back into ancestral environment at this stage on same axis
		}

		if withCue {
			for i, m := range indiv.Genome.E.Mat {
				for j, d := range m { //range over keys
					gproj += d * Gaxis.E.Mat[i][j]
				}
			}
		}
		if epig {
			for i, m := range indiv.Genome.F.Mat {
				for j, d := range m { //range over keys
					gproj += d * Gaxis.F.Mat[i][j]
				}
			}
		}
		for i, m := range indiv.Genome.G.Mat {
			for j, d := range m { //range over keys
				gproj += d * Gaxis.G.Mat[i][j]
			}
		}
		if hoc {
			for i, m := range indiv.Genome.Hg.Mat {
				for j, d := range m { //range over keys
					gproj += d * Gaxis.Hg.Mat[i][j]
				}
			}
			if hoi {
				for i, m := range indiv.Genome.Hh.Mat {
					for j, d := range m { //range over keys
						gproj += d * Gaxis.Hg.Mat[i][j]
					}
				}
			}
		}
		for i, m := range indiv.Genome.P.Mat {
			for j, d := range m { //range over keys
				gproj += d * Gaxis.P.Mat[i][j]
			}
		}
		/*
			for i, m := range indiv.Genome.Z.Mat {
				for j, d := range m { //range over keys
					gproj += d * Gaxis.Z.Mat[i][j]
				}
			}
		*/
		//fmt.Printf("Nov: %e\t Anc: %e\t G: %e\n", novpproj, ancpproj, gproj)
		fmt.Fprintf(fout, "%e\t %e\t %e\n", novpproj, ancpproj, gproj)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}
