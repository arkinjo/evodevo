package multicell

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"

	"gonum.org/v1/gonum/mat"
)

type Population struct { //Population of individuals
	Params  Settings
	Gen     int
	NovEnvs Cues //Novel Environment
	AncEnvs Cues // Ancestral Environment
	Indivs  []Indiv
}

type PopStats struct { // Statistics of population (mean values of quantities of interest)
	PEDot      float64 // (delta P, delta E) dot product
	PErr1      float64 // <|| p(e1) - e1 ||> (Nov)
	PErr0      float64 // <|| p(e0) - e0 ||> (Anc)
	PED10      float64 // <|| p(e1) - e0 ||> (Nov-Anc)
	PED01      float64 // <|| p(e0) - e1 ||> (Anc-Nov)
	WagFit     float64
	Fitness    float64
	Plasticity float64
	Div        float64
	NDevStep   float64
}

func (pop *Population) GetStats() PopStats {
	var stats PopStats
	mf := 0.0
	maxfit := 0.0
	merr1 := 0.0
	merr0 := 0.0
	md10 := 0.0
	md01 := 0.0
	ndev := 0
	mop := 0.0 // mean observed plasticity
	fn := float64(len(pop.Indivs))
	pa := NewCues(ncells, nenv)
	pv := NewCues(ncells, nenv)

	denv := 0.0
	for i, env := range pop.NovEnvs { //To normalize wrt change in environment cue
		denv += DistVecs1(env, pop.AncEnvs[i])
	}

	for _, indiv := range pop.Indivs {
		merr1 += indiv.Dp1e1
		merr0 += indiv.Dp0e0
		md10 += indiv.Dp1e0
		md01 += indiv.Dp0e1
		mf += indiv.Fit
		mop += indiv.Plasticity
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
	for i, p := range pa {
		for j := range p {
			pa[i][j] /= fn

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
	env0 := FlattenEnvs(pop.AncEnvs)
	env1 := FlattenEnvs(pop.NovEnvs)
	lenE := len(env1)
	dirE := NewVec(lenE)
	DiffVecs(dirE, env1, env0)
	NormalizeVec(dirE)

	mp1 := GetMeanVec(pop.GetFlatStateVec("P", 1))
	dirP := NewVec(lenE)
	DiffVecs(dirP, mp1, env0)
	NormalizeVec(dirP)

	stats.PEDot = DotVecs(dirP, dirE)
	stats.PErr1 = merr1 / fn
	stats.PErr0 = merr0 / fn
	stats.PED10 = md10 / fn
	stats.PED01 = md01 / fn
	meanfit := mf / fn
	stats.Fitness = meanfit
	stats.WagFit = meanfit / maxfit
	stats.NDevStep = float64(ndev) / fn
	stats.Plasticity = mop / (fn * denv)
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

	p := Population{default_settings, 0, envs0, envs1, indivs}
	return p
}

func (pop *Population) FromJSON(filename string) {
	pop.ClearGenome()
	fin, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}

	byteValue, _ := ioutil.ReadAll(fin)
	err = json.Unmarshal(byteValue, pop)
	if err != nil {
		log.Fatal(err)
	}

	err = fin.Close()
	if err != nil {
		log.Fatal(err)
	}
	log.Println("Successfully imported population from", filename)
}

func (pop *Population) ToJSON(filename string) {
	jsonpop, err := json.Marshal(pop) //JSON encoding of population as byte array
	if err != nil {
		log.Fatal(err)
	}
	fout, err := os.OpenFile(filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644) //create json file
	if err != nil {
		log.Fatal(err)
	}
	_, err = fout.Write(jsonpop)
	if err != nil {
		log.Fatal(err)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	log.Println("Successfuly exported population to", filename)
}

func (pop *Population) SetWagnerFitness() { //compute normalized fitness value similar to Wagner (1996).
	var mf float64
	for _, indiv := range pop.Indivs {
		if mf < indiv.Fit {
			mf = indiv.Fit
		}
	}
	for i, indiv := range pop.Indivs {
		pop.Indivs[i].WagFit = math.Max(indiv.Fit/mf, minWagnerFitness) //Zero fitness individuals that don't converge can still reproduce
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

	pop1 := NewPopulation(ncell, npop)
	pop1.Params = pop.Params
	pop1.Gen = pop.Gen
	pop1.NovEnvs = CopyCues(pop.NovEnvs)
	pop1.AncEnvs = CopyCues(pop.AncEnvs)
	for i, indiv := range pop.Indivs {
		pop1.Indivs[i] = indiv.Copy()
	}
	return pop1
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
		if withE {
			for i, m := range Gtilde.E.Mat {
				for j, v := range m {
					MeanGenome.E.Mat[i][j] += v * fnpop
				}
			}
		}
		if withF {
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
		if withH {
			for i, m := range Gtilde.H.Mat {
				for j, v := range m {
					MeanGenome.H.Mat[i][j] += v * fnpop
				}
			}
			if withJ {
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

func (pop *Population) Get_Environment_Axis() Cues { //CwithJce of axis defined using difference of environment cues
	axlength2 := 0.0

	e := pop.NovEnvs  //Cue in novel (present) environment
	e0 := pop.AncEnvs //Cue in ancestral (previous) environment
	v := NewVec(nenv + ncells)
	de := NewCues(ncells, nenv)

	for i, p := range e {
		DiffVecs(v, p, e0[i])
		axlength2 += Norm2Sq(v)
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

	new_population := Population{pop.Params, 0, pop.NovEnvs, pop.AncEnvs, nindivs} //resets embryonic values to zero!

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

	new_population := Population{pop.Params, 0, pop.NovEnvs, pop.AncEnvs, nindivs} //resets embryonic values to zero!

	return new_population
}
func (pop *Population) SortPopIndivs() {
	sort.Slice(pop.Indivs, func(i, j int) bool { return pop.Indivs[i].Id < pop.Indivs[j].Id })
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

	pop.SetWagnerFitness()
	//We might need a sorter here.
	pop.SortPopIndivs()

	return *pop
}

//Records population trajectory and writes files
func (pop0 *Population) Evolve(test bool, ftraj *os.File, jsonout string, nstep, epoch int) Population {
	pop := *pop0

	fmt.Fprintln(ftraj, "#Epoch\tGen\tNpop\tPhenoEnvDot \tMeanErr1 \tMeanErr0 \tMeanDp1e0 \tMeanDp0e1 \tFitness \tWag_Fit \tObs_Plas \tDiversity \tNdev") //header

	for istep := 1; istep <= nstep; istep++ {
		pop.DevPop(istep)
		if test {
			if jsonout != "" { //Export JSON population of each generation in test mode
				filename := fmt.Sprintf("%s_%3.3d.json", jsonout, pop.Gen)
				pop.ToJSON(filename)
			}
		}

		pstat := pop.GetStats()
		popsize := len(pop.Indivs)

		fmt.Fprintf(ftraj, "%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", epoch, istep, popsize, pstat.PEDot, pstat.PErr1, pstat.PErr0, pstat.PED10, pstat.PED01, pstat.Fitness, pstat.WagFit, pstat.Plasticity, pstat.Div, pstat.NDevStep)

		fmt.Printf("Evolve: %d\t<ME1>: %e\t<ME0>: %e\n", istep, pstat.PErr1, pstat.PErr0)

		pop = pop.PairReproduce(maxPop)
	}
	return pop
}

/*
func (pop *Population) Dump_Phenotypes(Filename string, gen int) {
	pop.DevPop(gen)

	trait := make([]float64, nenv) //trait part only

	fout, err := os.OpenFile(Filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
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
*/

func (pop *Population) Dump_Projections(Filename string, gen int, Gaxis Genome, Paxis Cues) {
	var ancpproj, novpproj, gproj float64
	pop.DevPop(gen) //Not needed for bugfixing

	anccphen := make(Vec, nenv+ncells)
	novcphen := make(Vec, nenv+ncells)

	mu := pop.Get_Mid_Env()
	Projfilename := fmt.Sprintf("%s_%3.3d.dat", Filename, gen)

	fout, err := os.OpenFile(Projfilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "#AncPhen \t NovPhen \t Genotype")

	for _, indiv := range pop.Indivs {
		ancpproj, novpproj, gproj = 0.0, 0.0, 0.0

		for i, env := range mu { //For each environment cue
			DiffVecs(anccphen, indiv.Bodies[IAncEnv].Cells[i].P, env) //centralize
			ancpproj += DotVecs(anccphen, Paxis[i])                   //Plot phenotype when pulled back into ancestral environment at this stage on same axis
		}

		for i, env := range mu { //For each environment cue
			DiffVecs(novcphen, indiv.Bodies[INovEnv].Cells[i].P, env) //centralize
			novpproj += DotVecs(novcphen, Paxis[i])
		}

		gproj = DotSpmats(indiv.Genome.G, Gaxis.G)
		gproj += DotSpmats(indiv.Genome.P, Gaxis.P)
		if withE {
			gproj += DotSpmats(indiv.Genome.E, Gaxis.E)
		}
		if withF {
			gproj += DotSpmats(indiv.Genome.F, Gaxis.F)
		}

		if withH {
			gproj += DotSpmats(indiv.Genome.H, Gaxis.H)
			if withJ {
				gproj += DotSpmats(indiv.Genome.J, Gaxis.J)
			}
		}

		fmt.Fprintf(fout, "%f\t %f\t %f\n", ancpproj, novpproj, gproj)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
}

func (pop *Population) GetFlatStateVec(istate string, ienv int) Dmat {
	vs0 := make([]Vec, 0)
	for _, indiv := range pop.Indivs {
		tv0 := make([]float64, 0)
		for _, cell := range indiv.Bodies[ienv].Cells {
			tv0 = append(tv0, cell.GetState(istate)...)
		}
		vs0 = append(vs0, tv0)
	}

	return vs0
}

func (pop *Population) GetFlatGenome() Dmat {
	vs := make([]Vec, 0)
	for _, indiv := range pop.Indivs {
		tv := indiv.Genome.FlatVec()
		vs = append(vs, tv)
	}
	return vs
}

func (pop *Population) GetPCA(state0 string, ienv0 int, state1 string, ienv1 int) (*mat.Dense, Vec, *mat.Dense) {
	s0 := pop.GetFlatStateVec(state0, ienv0)
	s1 := pop.GetFlatStateVec(state1, ienv1)

	_, _, ccmat := GetCrossCov(s0, s1, true, true)

	U, vals, V := GetSVD(ccmat)
	return U, vals, V
}
