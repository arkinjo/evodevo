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

	pcindexPtr := flag.Int("pcindex", -1, "index of principal component of environment")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 20, "magnitude of environmental change")
	tfilenamePtr := flag.String("traj_file", "traj.dat", "Filename of trajectories")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "Basenome of json output file")
	flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	multicell.SetSeedCue(int64(*seed_cuePtr))
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = *tfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.Omega = *omegaPtr

	ncells := *ncellsP
	pop0 := multicell.NewPopulation(ncells, *maxpopP) //with randomized genome to start
	maxpcindex := ncells*(multicell.GetNenv()+ncells) - 1
	pcindex := multicell.MinInt(*pcindexPtr, maxpcindex)

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

	if pcindex < 0 { //If no directions given
		NovEnvs = multicell.ChangeEnvs(OldEnvs, denv) //Randomize
	} else { //If directions are given
		NovEnvs = multicell.PCA2Cue(&pop0, pcindex)
		fmt.Println("NovEnv(PCA)", pcindex, ":", NovEnvs)
	}

	popstart.NovEnvs = NovEnvs //control size of perturbation of environment cue vector at start of epoch.

	tevol := time.Now()

	fmt.Println("Evolving in novel environment :", popstart.NovEnvs)
	fmt.Println("Ancestral environment :", popstart.AncEnvs)

	pop1 := popstart.Evolve(true, ftraj, json_out, epochlength, 1)
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
