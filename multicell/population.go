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
			//fmt.Println("Id:", indiv.Id, "Fitness:", indiv.Fit, "Current Maximum fitness value:", mf)
		}
	}
	for i, indiv := range pop.Indivs {
		pop.Indivs[i].WagFit = math.Max(indiv.Fit/mf, minWagnerFitness) //Zero fitness individuals that don't converge can still reproduce
		//fmt.Println("Id:", indiv.Id, "Fitness:", indiv.Fit, "Wagner Fitness:", indiv.WagFit)
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
	ncell := len(pop.Indivs[0].Copy().Copies[INovEnv].Ctypes) //number of cells
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

func (pop *Population) GetMeanFitness() (float64, float64) { //average fitness of population
	mf := 0.0
	maxfit := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mf += indiv.Fit
		if indiv.Fit > maxfit {
			maxfit = indiv.Fit
		}
	}
	meanfit := mf / fn
	wagfit := meanfit / maxfit

	return meanfit, wagfit
}

/*
func (pop *Population) GetMeanWagFit() float64 { //average Wagner fitness of population
	mf := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mf += indiv.WagFit
	}

	return mf / fn
}
*/

func (pop *Population) GetMSE() float64 { //average observed plasticity of population
	mse := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mse += indiv.MSE
	}

	return mse / fn
}

func (pop *Population) GetMeanObsPlasticity() float64 { //average observed plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	denv := 0.0
	for i, env := range pop.Envs { //To normalize wrt change in environment cue
		denv += Hammingdist(env, pop.RefEnvs[i])
	}
	for _, indiv := range pop.Indivs {
		mp += indiv.ObsPlas
	}

	return mp / (fn * denv)
}

func (pop *Population) GetMeanAncCuePlasticity() float64 { //average cue plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		//fmt.Println("ID:", indiv.Id, "CuePlas:", indiv.CuePlas)
		mp += indiv.AncCuePlas
	}

	return mp / fn
}

func (pop *Population) GetMeanNovCuePlasticity() float64 { //average cue plasticity of population
	mp := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		//fmt.Println("ID:", indiv.Id, "CuePlas:", indiv.CuePlas)
		mp += indiv.NovCuePlas
	}

	return mp / fn
}

/*
func (pop *Population) GetMeanUtility() float64 { //average utility of population
	mu := 0.0
	fn := float64(len(pop.Indivs))
	for _, indiv := range pop.Indivs {
		mu += indiv.Util
	}

	return mu / fn
}
*/

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
		for _, cell := range ind.Copies[INovEnv].Ctypes {
			cv = append(cv, cell.P)
		}
	}
	div := GetCueVar(cv)

	return div
}

/*
func (pop *Population) GetGeneFunc() float64 {
	gf := 0.0
	for _, indiv := range pop.Indivs {
		gf += indiv.MSE
		gf += indiv.ObsPlas
	}
	return gf
}
*/

func (pop *Population) GetMeanPhenotype(gen int) Cues { //elementwise average phenotype of population; output as slice instead of cue struct
	npop := len(pop.Indivs)
	MeanPhenotype := NewCues(ncells, nenv)
	pop.DevPop(gen)

	for _, indiv := range pop.Indivs {
		for i, c := range indiv.Copies[INovEnv].Ctypes {
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

func (pop *Population) Selection(nNewPop int) []Indiv { //Selects parents for new population
	npop := len(pop.Indivs)
	//var parents []Indiv //Does this even work?
	parents := make([]Indiv, 0)
	ipop := 0
	cnt := 0
	for ipop < nNewPop && cnt < 1000*nNewPop {
		cnt += 1
		k := rand.Intn(npop)
		ind := pop.Indivs[k]
		r := rand.Float64()
		if r < ind.WagFit {
			parents = append(parents, ind)
			ipop += 1
		}
	}
	return parents
}

func (pop *Population) Reproduce(nNewPop int) Population { //Crossover
	parents := pop.Selection(nNewPop)
	nindivs := make([]Indiv, 0)
	npop := len(parents)

	//fmt.Println("Number of individuals that survived selection :", npop)

	//ipop := 0
	//cnt := 0
	for len(nindivs) < nNewPop { //Randomly reproduce among survivors
		k := rand.Intn(npop)
		l := rand.Intn(npop)
		dad := parents[k]
		mom := parents[l]
		kid0, kid1 := Mate(&dad, &mom)
		nindivs = append(nindivs, kid0)
		nindivs = append(nindivs, kid1)
		//ipop += 2
	}

	for i := range nindivs {
		nindivs[i].Id = i //Relabels individuals according to position in array
	}

	new_population := Population{0, pop.Envs, pop.RefEnvs, nindivs} //resets embryonic values to zero!

	return new_population

}

func (pop *Population) PairReproduce(nNewPop int) Population { //Crossover in ordered pairs; as in Wagner's
	var index int
	parents := pop.Selection(nNewPop)
	//nparents := len(parents)
	nindivs := make([]Indiv, 0)

	for index < nNewPop && len(nindivs) < nNewPop { //Forced reproduction in ordered pairs; may cause bugs when population has an odd number of survivors
		dad := parents[index]
		mom := parents[index+1]
		kid0, kid1 := Mate(&dad, &mom)
		nindivs = append(nindivs, kid0)
		nindivs = append(nindivs, kid1)
		index = len(nindivs) //update
	}

	for i := range nindivs {
		nindivs[i].Id = i //Relabels individuals according to position in array
	}

	new_population := Population{0, pop.Envs, pop.RefEnvs, nindivs} //resets embryonic values to zero!

	return new_population
}

/*
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
*/

func (pop *Population) DevPop(gen int) Population {
	novenv := NewCues(ncells, nenv)
	ancenv := NewCues(ncells, nenv)

	pop.Gen = gen

	ch := make(chan Indiv) //channels for parallelization
	for _, indiv := range pop.Indivs {
		go func(indiv Indiv) {
			copy(novenv, pop.Envs)
			copy(ancenv, pop.RefEnvs)
			ch <- indiv.CompareDev(novenv, ancenv)
		}(indiv)
	}
	for i := range pop.Indivs {
		pop.Indivs[i] = <-ch //Update output results
	}

	/*
		for i, indiv := range pop.Indivs {
			copy(novenv, pop.Envs)
			copy(ancenv, pop.RefEnvs)
			pop.Indivs[i] = indiv.CompareDev(novenv, ancenv)
		}
	*/

	//We might need a sorter here.
	pop.SortPopIndivs()

	return *pop
}

func Evolve(test bool, tfilename, jsonout, gidfilename string, nstep, epoch int, init_pop *Population) Population { //Records population trajectory and writes files
	var jfilename, id_filename, id, dadid, momid string
	//var Fitness, CuePlas, ObsPlas, Polyp, Div, Util float64
	var MSE, WagFit, Fitness, AncCuePlas, NovCuePlas, ObsPlas, Polyp, Div float64
	var popsize int

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

	//track := true
	//MSE0 := 0.0
	//EMA_MSE := 0.0

	//Pl0 := 0.0
	//dPl := 0.0
	//EMA_Pl := 0.0

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

		pop.SetWagnerFitness()
		/*
			for _, indiv := range pop.Indivs {
				fmt.Println("Before reproduction: Id: ", indiv.Id, "Wagner Fitness: ", indiv.WagFit, "Fitness: ", indiv.Fit)
			}
		*/

		Fitness, WagFit = pop.GetMeanFitness()
		MSE = pop.GetMSE()
		AncCuePlas = pop.GetMeanAncCuePlasticity()
		NovCuePlas = pop.GetMeanNovCuePlasticity()
		ObsPlas = pop.GetMeanObsPlasticity()
		Polyp = pop.GetMeanPp()
		Div = pop.GetDiversity()
		//Util = pop.GetMeanUtility()
		popsize = len(pop.Indivs)
		/*
			EMA_MSE = 2.0/(l_EMA+1.0)*MSE + (1-2/(l_EMA+1.0))*EMA_MSE
			EMA_Pl = 2.0/(l_EMA+1.0)*Pl + (1-2/(l_EMA+1.0))*EMA_Pl
		*/
		fout, err := os.OpenFile(tfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}

		//fmt.Fprintf(fout, "%d\t%d\t%f\t%e\t%e\t%e\t%e\t%e\n", epoch, istep, Fitness, CuePlas, ObsPlas, Polyp, Div, Util)

		//fmt.Fprintf(fout, "%d\t%d\t%f\t%e\t%e\t%e\t%e\n", epoch, istep, MSE, Fitness, Pl, Polyp, Div, Util)
		fmt.Fprintf(fout, "%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", epoch, istep, popsize, MSE, Fitness, WagFit, AncCuePlas, NovCuePlas, ObsPlas, Polyp, Div)

		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

		//fmt.Printf("Evol_step: %d\t <Fit>: %f\t <Pl>:%e\t <Pp>:%e\t <Div>:%e \t <u>:%e\n ", istep, Fitness, Pl, Polyp, Div, Util)
		fmt.Printf("Evol_step: %d\t <Npop>: %d\t <MSE>: %e\t <Fit>: %e\t <WFit>: %e\t <ACPl>:%e\t <NCPl>:%e\t <OPl>:%e\t <Pp>:%e\t <Div>:%e \n ", istep, popsize, MSE, Fitness, WagFit, AncCuePlas, NovCuePlas, ObsPlas, Polyp, Div)

		pop = pop.PairReproduce(maxPop)

		/*
			for _, indiv := range pop.Indivs {
				fmt.Println("After reproduction: Id: ", indiv.Id, "Wagner Fitness: ", indiv.WagFit, "Fitness: ", indiv.Fit) //Bugtest
			}
		*/

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

func (pop *Population) Dump_Phenotypes(Filename string, gen int) {
	pop.DevPop(gen)

	trait := make([]float64, nenv) //trait part only

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	for _, indiv := range pop.Indivs {
		for _, cell := range indiv.Copies[INovEnv].Ctypes {
			copy(trait, cell.P)
			for _, v := range trait {
				fmt.Fprintf(fout, "%f\t", v)
			}
		}
		fmt.Fprint(fout, "\n") //Each line is phenotype vector concatenation of one individual
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

func (pop *Population) Dump_Projections(Filename string, gen int, Gaxis Genome, Paxis Cues) {
	var defpproj, ancpproj, novpproj, gproj float64
	pop.DevPop(gen) //Not needed for bugfixing

	anccphen := make(Vec, nenv+ncells)
	novcphen := make(Vec, nenv+ncells)
	defcphen := make(Vec, nenv+ncells)

	mu := pop.Get_Mid_Env()
	//fmt.Println("Middle environment : ", mu)
	Projfilename := fmt.Sprintf("../analysis/%s_%d.dat", Filename, gen)

	fout, err := os.OpenFile(Projfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "DefPhen \t AncPhen \t NovPhen \t Genotype")

	for _, indiv := range pop.Indivs {
		defpproj, ancpproj, novpproj, gproj = 0.0, 0.0, 0.0, 0.0

		for i, env := range mu { //For each environment cue
			diffVecs(defcphen, indiv.Copies[INoEnv].Ctypes[i].P, env) //centralize
			defpproj += Innerproduct(defcphen, Paxis[i])
		}

		for i, env := range mu { //For each environment cue
			diffVecs(anccphen, indiv.Copies[IAncEnv].Ctypes[i].P, env) //centralize
			ancpproj += Innerproduct(anccphen, Paxis[i])               //Plot phenotype when pulled back into ancestral environment at this stage on same axis
		}

		for i, env := range mu { //For each environment cue
			diffVecs(novcphen, indiv.Copies[INovEnv].Ctypes[i].P, env) //centralize
			novpproj += Innerproduct(novcphen, Paxis[i])
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
		fmt.Fprintf(fout, "%f\t %f\t %f\t %f\n", defpproj, novpproj, ancpproj, gproj)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}
