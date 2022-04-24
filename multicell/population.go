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
	NovEnvs Cues //Novel Environment
	AncEnvs Cues // Ancestral Environment
	Indivs  []Indiv
}

type PopStats struct { // Statistics of population (mean values of quantities of interest)
	CDist      float64 // Distance between average phenotype and env.
	MeanErr1   float64 // || p(e1) - e1 || (Nov)
	MeanErr0   float64 // || p(e0) - e0 || (Anc)
	MeanDiff   float64 // || p(e1) - e0 || (Nov-Anc)
	WagFit     float64
	Fitness    float64
	AncCuePlas float64
	NovCuePlas float64
	ObsPlas    float64
	Polyp      float64
	Div        float64
	NDevStep   float64
}

func (pop *Population) GetStats() PopStats {
	var stats PopStats
	mf := 0.0
	maxfit := 0.0
	merr1 := 0.0
	merr0 := 0.0
	mdiff := 0.0
	ndev := 0
	mop := 0.0 // mean observed plasticity
	ap := 0.0  // mean ancestral plasticity
	np := 0.0  // mean (novel) plasticity
	pp := 0.0  // polyphenism
	fn := float64(len(pop.Indivs))
	pa := NewCues(ncells, nenv)
	pv := NewCues(ncells, nenv)

	denv := 0.0
	for i, env := range pop.NovEnvs { //To normalize wrt change in environment cue
		denv += Hammingdist(env, pop.AncEnvs[i])
	}

	for _, indiv := range pop.Indivs {
		merr1 += indiv.getMeanErr(INovEnv)
		merr0 += indiv.getMeanErr(IAncEnv)
		mdiff += indiv.Dp1e0
		mf += indiv.Fit
		mop += indiv.ObsPlas
		ap += indiv.AncCuePlas
		np += indiv.NovCuePlas
		pp += indiv.Pp
		ndev += indiv.Bodies[INovEnv].NDevStep

		if indiv.Fit > maxfit {
			maxfit = indiv.Fit
		}
		for i, cell := range indiv.Bodies[INovEnv].Cells {
			for j, t := range cell.P {
				pa[i][j] += t
			}
		}

	}
	cdist := 0.0
	for i, p := range pa {
		for j := range p {
			pa[i][j] /= fn
			d := pop.NovEnvs[i][j] - pa[i][j]
			cdist += d * d
		}
	}
	for _, indiv := range pop.Indivs {
		for i, cell := range indiv.Bodies[INovEnv].Cells {
			for j, t := range cell.P {
				d := pa[i][j] - t
				pv[i][j] += d * d
			}
		}
	}

	div := 0.0
	for _, pi := range pv {
		for _, t := range pi {
			div += t / fn
		}
	}

	stats.CDist = math.Sqrt(cdist)
	stats.MeanErr1 = merr1 / fn
	stats.MeanErr0 = merr0 / fn
	stats.MeanDiff = mdiff / fn
	meanfit := mf / fn
	stats.Fitness = meanfit
	stats.WagFit = meanfit / maxfit
	stats.NDevStep = float64(ndev) / fn
	stats.ObsPlas = mop / (fn * denv)
	stats.AncCuePlas = ap / fn
	stats.NovCuePlas = np / fn
	stats.Polyp = pp / fn
	stats.Div = div

	return stats
}

func NewPopulation(ncell, npop int) Population { //Initialize new population
	envs0 := NewCues(ncell, nenv)
	envs1 := NewCues(ncell, nenv)

	indivs := make([]Indiv, npop)
	for i := range indivs {
		indivs[i] = NewIndiv(i)
	}

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
	ncell := len(pop.Indivs[0].Bodies[INovEnv].Cells) //number of cells
	//fmt.Println("Copying ",ncell,"-cell individuals")
	pop1 := NewPopulation(ncell, npop)
	pop1.Gen = pop.Gen
	pop1.NovEnvs = CopyCues(pop.NovEnvs)
	pop1.AncEnvs = CopyCues(pop.AncEnvs)
	for i, indiv := range pop.Indivs {
		pop1.Indivs[i] = indiv.Copy()
	}
	return pop1
}

func (pop *Population) SortPopIndivs() {
	sort.Slice(pop.Indivs, func(i, j int) bool { return pop.Indivs[i].Id < pop.Indivs[j].Id }) //Hopefully this works
}

func (pop *Population) GetMeanPhenotype(gen int) Cues { //elementwise average phenotype of population; output as slice instead of cue struct
	npop := len(pop.Indivs)
	MeanPhenotype := NewCues(ncells, nenv)
	pop.DevPop(gen)

	for _, indiv := range pop.Indivs {
		for i, c := range indiv.Bodies[INovEnv].Cells {
			for j, p := range c.P {
				MeanPhenotype[i][j] += p / float64(npop)
			}
		}
	}
	return MeanPhenotype
}

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
			for i, m := range Gtilde.H.Mat {
				for j, v := range m {
					MeanGenome.H.Mat[i][j] += v * fnpop
				}
			}
			if hoi {
				for i, m := range Gtilde.J.Mat {
					for j, v := range m {
						MeanGenome.J.Mat[i][j] += v * fnpop
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

	e := pop.NovEnvs  //Cue in novel (present) environment
	e0 := pop.AncEnvs //Cue in ancestral (previous) environment
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
	e := pop.NovEnvs  // novel environment
	e0 := pop.AncEnvs // ancestral environment

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

	new_population := Population{0, pop.NovEnvs, pop.AncEnvs, nindivs} //resets embryonic values to zero!

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

	new_population := Population{0, pop.NovEnvs, pop.AncEnvs, nindivs} //resets embryonic values to zero!

	return new_population
}

func (pop *Population) DevPop(gen int) Population {
	pop.Gen = gen

	ch := make(chan Indiv) //channels for parallelization
	for _, indiv := range pop.Indivs {
		go func(indiv Indiv) {
			ch <- indiv.Develop(pop.AncEnvs, pop.NovEnvs)
		}(indiv)
	}
	for i := range pop.Indivs {
		pop.Indivs[i] = <-ch //Update output results
	}

	//We might need a sorter here.
	pop.SortPopIndivs()

	return *pop
}

func Evolve(test bool, tfilename, jsonout, gidfilename string, nstep, epoch int, init_pop *Population) Population { //Records population trajectory and writes files
	var jfilename, id_filename, id, dadid, momid string
	//var Fitness, CuePlas, ObsPlas, Polyp, Div, Util float64
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

		pstat := pop.GetStats()
		popsize = len(pop.Indivs)

		fout, err := os.OpenFile(tfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}

		fmt.Fprintf(fout, "%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", epoch, istep, popsize, pstat.CDist, pstat.MeanErr1, pstat.MeanErr0, pstat.MeanDiff, pstat.Fitness, pstat.WagFit, pstat.AncCuePlas, pstat.NovCuePlas, pstat.ObsPlas, pstat.Polyp, pstat.Div, pstat.NDevStep)

		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf("Evol_step: %d\t <Npop>: %d\t<CD>: %e\t<ME1>: %e\t<ME0>: %e\t<MDf>: %e\t<Fit>: %e\t<WFit>: %e\t<ACPl>: %e\t<NCPl>: %e\t<OPl>: %e\t<Pp>: %e\t<Div>: %e\t<Ndev>: %e\n ", istep, popsize, pstat.CDist, pstat.MeanErr1, pstat.MeanErr0, pstat.MeanDiff, pstat.Fitness, pstat.WagFit, pstat.AncCuePlas, pstat.NovCuePlas, pstat.ObsPlas, pstat.Polyp, pstat.Div, pstat.NDevStep)

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
		for _, cell := range indiv.Bodies[INovEnv].Cells {
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
			diffVecs(defcphen, indiv.Bodies[INoEnv].Cells[i].P, env) //centralize
			defpproj += Innerproduct(defcphen, Paxis[i])
		}

		for i, env := range mu { //For each environment cue
			diffVecs(anccphen, indiv.Bodies[IAncEnv].Cells[i].P, env) //centralize
			ancpproj += Innerproduct(anccphen, Paxis[i])              //Plot phenotype when pulled back into ancestral environment at this stage on same axis
		}

		for i, env := range mu { //For each environment cue
			diffVecs(novcphen, indiv.Bodies[INovEnv].Cells[i].P, env) //centralize
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
			for i, m := range indiv.Genome.H.Mat {
				for j, d := range m { //range over keys
					gproj += d * Gaxis.H.Mat[i][j]
				}
			}
			if hoi {
				for i, m := range indiv.Genome.J.Mat {
					for j, d := range m { //range over keys
						gproj += d * Gaxis.H.Mat[i][j]
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
		//gproj = gproj / math.Sqrt(float64(genelength*ngenes)) //Normalize wrt genome size
		fmt.Fprintf(fout, "%f\t %f\t %f\t %f\n", defpproj, ancpproj, novpproj, gproj)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}
