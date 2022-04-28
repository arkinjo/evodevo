package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var PG_Filename string  //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var PCA_Filename string
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
	seed_cuePtr := flag.Int("seed_cue", 1, "random seed for environmental cue")
	//epochPtr := flag.Int("nepoch", 1, "number of epochs")
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case

	pcindexPtr := flag.Int("pcindex", 0, "index of principal component of environment")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 20, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename", "traj", "name of file of trajectories")
	pcafilenamePtr := flag.String("pcafilename", "", "name of file of principal trait vector concatanation")
	pgfilenamePtr := flag.String("pgfilename", "", "name of file of projected phenotypes and genotypes") //default to empty string
	gidfilenamePtr := flag.String("gidfilename", "", "name of file of geneology of ids")                 //default to empty string
	jsoninPtr := flag.String("jsonin", "", "json file of input population")                              //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")
	//testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
	flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	multicell.SetSeedCue(int64(*seed_cuePtr))
	//maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = fmt.Sprintf("../analysis/%s.dat", *tfilenamePtr)
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr
	PCA_Filename = *pcafilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.Omega = *omegaPtr

	pcindex := multicell.MinInt(*pcindexPtr, multicell.GetNcells()*multicell.GetNenv()-1)

	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start

	if json_in != "" { //read input population as a json file, if given
		pop0.FromJSON(json_in)
		multicell.SetParams(pop0.Params)
	} else {
		log.Fatal("JSON input required.")
	}

	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	ftraj, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0.Copy()
	AncEnvs := multicell.CopyCues(pop0.NovEnvs)
	OldEnvs := multicell.CopyCues(pop0.NovEnvs)
	popstart.AncEnvs = AncEnvs
	NovEnvs := multicell.NewCues(multicell.GetNcells(), multicell.GetNenv()) //Declaration
	if PCA_Filename == "" {                                                  //If no directions given
		NovEnvs = multicell.ChangeEnvs(OldEnvs, denv) //Randomize
	} else { //If directions are given
		//pcafilename := fmt.Sprintf("%s.dat", PCA_Filename)
		pcacues, pcavecs := multicell.PCAtoCue(PCA_Filename)
		NovEnvs = multicell.CopyCues(pcacues[pcindex]) //copy
		fmt.Println(pcavecs)
	}

	popstart.NovEnvs = NovEnvs //control size of perturbation of environment cue vector at start of epoch.

	tevol := time.Now()

	fmt.Println("Evolving in novel environment :", popstart.NovEnvs)
	fmt.Println("Ancestral environment :", popstart.AncEnvs)
	gidfilename := fmt.Sprintf("%s_full", Gid_Filename)

	pop1 := popstart.Evolve(true, ftraj, json_out, gidfilename, epochlength, 1)
	fmt.Println("End of epoch")

	dtevol := time.Since(tevol)
	fmt.Println("Time taken to simulate evolution :", dtevol)
	fmt.Println(pop1.NovEnvs)
	err = ftraj.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Trajectory of population written to", T_Filename)
	fmt.Printf("JSON encoding of populations written to %s_*.json \n", json_out)

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
