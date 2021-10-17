package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var PG_Filename string  //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"

var CopyequalGs float64 //bugtesting variable
var JSONequalGs float64 //bugtesting variable
var DevequalGs float64  //bugtesting variable

//var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	seedPtr := flag.Int("seed", 1, "random seed")
	//epochPtr := flag.Int("nepoch", 1, "number of epochs")
	maxpopsizePtr := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuestrengthPtr := flag.Float64("cuestrength", 1.0, "control size of var contribution of environmental cue")
	hoistrengthPtr := flag.Float64("hoistrength", 1.0, "control size of var contribution of higher order interactions")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	//HOIPtr := flag.Bool("HOI", true, "Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 2, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename", "traj", "name of file of trajectories")
	pgfilenamePtr := flag.String("pgfilename", "", "name of file of projected phenotypes and genotypes") //default to empty string
	gidfilenamePtr := flag.String("gidfilename", "", "name of file of geneology of ids")                 //default to empty string
	jsoninPtr := flag.String("jsonin", "", "json file of input population")                              //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")
	//testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
	flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	//maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = fmt.Sprintf("../analysis/%s.dat", *tfilenamePtr)
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.Omega = *omegaPtr

	multicell.SetMaxPop(*maxpopsizePtr)
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *epigPtr, *HOCPtr)

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop()) //with randomized genome to start

	if json_in != "" { //read input population as a json file, if given
		fmt.Println("Importing initial population")
		jfilename := fmt.Sprintf("../pops/%s.json", json_in)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}

		byteValue, _ := ioutil.ReadAll(popin)
		pop0.ClearGenome() //Clear genome before unmarshalling json
		err = json.Unmarshal(byteValue, &pop0)
		if err != nil {
			log.Fatal(err)
		}

		err = popin.Close()
		if err != nil {
			log.Fatal(err)
		}
		fmt.Println("Successfully imported population")
	} else {
		fmt.Println("Randomizing initial population")
		pop0.RandomizeGenome()
	}

	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	fout, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}

	fmt.Fprintln(fout, "Epoch \t Generation \t Fitness \t Cue_Plas \t Obs_Plas \t Polyphenism \t Diversity \t Utility") //header
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0.Copy()
	AncEnvs := multicell.CopyCues(pop0.Envs)
	OldEnvs := multicell.CopyCues(pop0.Envs)
	popstart.RefEnvs = AncEnvs
	NovEnvs := multicell.ChangeEnvs(OldEnvs, denv)
	popstart.Envs = NovEnvs //control size of perturbation of environment cue vector at start of epoch.

	tevol := time.Now()

	fmt.Println("Evolving in novel environment :", popstart.Envs)
	fmt.Println("Ancestral environment :", popstart.RefEnvs)
	gidfilename := fmt.Sprintf("%s_full", Gid_Filename)

	pop1 := multicell.Evolve(true, T_Filename, json_out, gidfilename, epochlength, 1, &popstart)
	fmt.Println("End of epoch")

	dtevol := time.Since(tevol)
	fmt.Println("Time taken to simulate evolution :", dtevol)

	fmt.Println("Dumping projections")
	tdump := time.Now()
	pop := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop())
	g0 := pop0.GetMeanGenome()
	g1 := pop1.GetMeanGenome()
	//dG := multicell.TestEqualGenomes(g0, g1)
	//fmt.Println("Difference in genome:", dG)

	Gaxis := multicell.NewGenome()
	//ZeroGenome := multicell.NewGenome()
	multicell.DiffGenomes(&Gaxis, &g1, &g0) //Pointers for bugfix; don't ask why; it just works!
	//dG = multicell.TestEqualGenomes(ZeroGenome, Gaxis)
	//fmt.Println("Axis length:", dG) //Should be equal to difference in genome; currently this is stuck at zero value of genome.

	Gaxis = Gaxis.NormalizeGenome()
	//dG = multicell.TestEqualGenomes(ZeroGenome, Gaxis)
	//fmt.Println("Normalized Axis length:", dG)

	Paxis := pop1.Get_Environment_Axis() //Measure everything in direction of ancestral -> novel environment
	//fmt.Println("Mean environment:", pop1.Get_Mid_Env())
	//fmt.Println("Change in environment proportional to:", Paxis)

	//bugfixpop := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop) //Bugfixing

	for gen := 1; gen <= epochlength; gen++ { //Also project population after pulling back to ancestral environment.
		jfilename := fmt.Sprintf("../pops/%s_%d.json", json_out, gen)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}

		pop.ClearGenome()
		byteValue, _ := ioutil.ReadAll(popin)
		err = json.Unmarshal(byteValue, &pop)
		if err != nil {
			log.Fatal(err)
		}

		//fmt.Println("Mean environment:", pop.Get_Mid_Env())
		//fmt.Println("Change in environment proportional to:", pop.Get_Environment_Axis())
		//dG = multicell.TestEqualPopGenomes(bugfixpop, pop)
		//fmt.Println("Norm of imported genotype :", dG)

		err = popin.Close()
		if err != nil {
			log.Fatal(err)
		}

		pop.Dump_Projections(PG_Filename, gen, Gaxis, Paxis)
	}
	dtdump := time.Since(tdump)
	fmt.Println("Time taken to dump projections :", dtdump)
	fmt.Println("Making DOT genealogy file")
	tdot := time.Now()
	nanctraj := multicell.DOT_Genealogy(Gid_Filename, json_out, epochlength, multicell.GetMaxPop())
	//fmt.Println(nanctraj)
	dtdot := time.Since(tdot)
	fmt.Println("Time taken to make dot file :", dtdot)
	fmt.Println("Dumping number of ancestors")
	nancfilename = fmt.Sprintf("../analysis/%s_nanc.dat", Gid_Filename)
	fout, err = os.OpenFile(nancfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "Generation \t Ancestors")
	for i, n := range nanctraj {
		fmt.Fprintf(fout, "%d\t%d\n", i+1, n)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Trajectory of population written to", T_Filename)
	fmt.Printf("Projections written to %s_*.dat \n", PG_Filename)
	fmt.Printf("Genealogy of final generation written to %s.dot\n", Gid_Filename)
	fmt.Printf("Number of ancestors of final generation written to %s\n", nancfilename)
	fmt.Printf("JSON encoding of populations written to %s_*.json \n", json_out)

	/*
		//fmt.Println("Trajectory of environment :", envtraj)

		//if !CopyequalGs {
		fmt.Println("Error in copying genomes :", CopyequalGs)
		//}
		//if !DevequalGs {
		fmt.Println("Genome error in development:", DevequalGs) //This is due to misordering after development; individuals are not developed in same order.
		//}
		//if !JSONequalGs {
		fmt.Println("Genome error in json:", JSONequalGs)
		//}
	*/

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
