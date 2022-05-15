package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj.dat"
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"
var jfilename string

func main() {
	t0 := time.Now()
	testP := flag.Bool("test", false, "Test run or not")
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	elayerP := flag.Bool("layerE", true, "Environmental cue layer")
	flayerP := flag.Bool("layerF", true, "Epigenetic layer")
	hlayerP := flag.Bool("layerH", true, "Higher order complexes")
	jlayerP := flag.Bool("layerJ", true, "Interactions in higher order interactions")
	pfbackP := flag.Bool("pfback", true, "Phenotype feedback to input")
	ncellsP := flag.Int("ncells", 1, "Number of cell types")
	sdNoiseP := flag.Float64("sdNoise", 0.05, "Std.Dev. of environmental noise")

	seedPtr := flag.Int("seed", 13, "random seed")
	seed_cuePtr := flag.Int("seed_cue", 7, "random seed for environmental cue")
	epochPtr := flag.Int("nepoch", 20, "number of epochs")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")

	denvPtr := flag.Int("denv", 20, "magnitude of environmental change")
	tfilenamePtr := flag.String("traj_file", "traj.dat", "filename of trajectories")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")
	flag.Parse()

	var settings = multicell.Settings{*maxpopP, *ncellsP, *elayerP, *flayerP, *hlayerP, *jlayerP, *pfbackP, *sdNoiseP}

	log.Println("seed=", *seedPtr, "seed_cue=", *seed_cuePtr)
	multicell.SetSeed(int64(*seedPtr))
	multicell.SetSeedCue(int64(*seed_cuePtr))

	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = *tfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	test_flag := *testP

	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP)
	pop0.Params = settings

	if json_in != "" { //read input population as a json file, if given
		pop0.FromJSON(json_in)
	} else {
		fmt.Println("Randomizing initial population")
		pop0.RandomizeGenome()
	}
	multicell.SetParams(pop0.Params)

	ftraj, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0
	popstart.NovEnvs = multicell.RandomEnvs(multicell.GetNcells(), multicell.GetNenv(), 0.5)
	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	envtraj := make([]multicell.Cues, 1) //Trajectory of environment cue
	envtraj[0] = popstart.AncEnvs
	novvec := make([]bool, 0)

	log.Println("AncEnvs", 0, ":", popstart.AncEnvs)
	for epoch := 1; epoch <= maxepochs; epoch++ {
		tevol := time.Now()
		log.Println("NovEnvs", epoch, ":", popstart.NovEnvs)
		envtraj = append(envtraj, popstart.NovEnvs)
		if epoch != 0 {
			fmt.Println("Epoch ", epoch, "has environments", popstart.NovEnvs)
		}

		pop1 := popstart.Evolve(test_flag, ftraj, json_out, epochlength, epoch)
		fmt.Println("End of epoch", epoch)

		if !test_flag && epoch == maxepochs { //Export output population; just before epoch change
			//Update to environment just before epoch change
			pop1.AncEnvs = multicell.CopyCues(pop1.NovEnvs)
			pop1.ToJSON(json_out)
		}
		dtevol := time.Since(tevol)
		fmt.Println("Time taken to simulate evolution :", dtevol)

		popstart = pop1 //Update population after evolution.

		OldEnvs := multicell.CopyCues(popstart.NovEnvs)
		popstart.AncEnvs = OldEnvs
		popstart.NovEnvs = multicell.ChangeEnvs2(OldEnvs, denv)
		err = multicell.DeepVec3NovTest(popstart.NovEnvs, envtraj)
		if err != nil {
			fmt.Println(err)
		}
		novvec = append(novvec, err == nil)
	}
	err = ftraj.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Trajectory of population written to %s \n", T_Filename)
	fmt.Printf("JSON encoding of evolved population written to %s \n", jfilename)

	fmt.Println("Novelty of environment cue :", novvec)
	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
