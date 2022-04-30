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
var P_Filename string = "pvec"
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"
var jfilename string

func main() {
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")

	elayerP := flag.Bool("elayer", true, "Environmental cue layer")
	flayerP := flag.Bool("flayer", true, "Epigenetic layer")
	hlayerP := flag.Bool("hlayer", true, "Higher order complexes")
	jlayerP := flag.Bool("jlayer", true, "Interactions in higher order interactions")
	pfbackP := flag.Bool("pfback", true, "Phenotype feedback to input")
	ncellsP := flag.Int("ncells", 1, "Number of cell types")
	sdNoiseP := flag.Float64("sdNoise", 0.05, "Std.Dev. of environmental noise")

	seedPtr := flag.Int("seed", 1, "random seed")
	seed_cuePtr := flag.Int("seed_cue", 1, "random seed for environmental cue")
	epochPtr := flag.Int("nepoch", 20, "number of epochs")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")

	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 20, "magnitude of environmental change")
	tfilenamePtr := flag.String("traj_file", "traj", "filename of trajectories")
	pfilenamePtr := flag.String("pheno_file", "pvec", "filename of phenotypes")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")
	flag.Parse()

	var settings = multicell.Settings{*maxpopP, *ncellsP, *elayerP, *flayerP, *hlayerP, *jlayerP, *pfbackP, *sdNoiseP}

	multicell.SetSeed(int64(*seedPtr))
	multicell.SetSeedCue(int64(*seed_cuePtr))

	multicell.SetParams(settings)

	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = *tfilenamePtr
	P_Filename = fmt.Sprintf("../analysis/%s.dat", *pfilenamePtr) //all phenotypes go to analysis directory
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.Omega = *omegaPtr

	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP)
	pop0.Params = settings

	if json_in != "" { //read input population as a json file, if given
		pop0.FromJSON(json_in)
	} else {
		fmt.Println("Randomizing initial population")
		pop0.RandomizeGenome()
	}

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

	for epoch := 1; epoch <= maxepochs; epoch++ {
		tevol := time.Now()
		envtraj = append(envtraj, popstart.NovEnvs)
		//existing envtraj entries should not be updated with each append/update.
		//Could it be reading popstart.Envs on each append?
		//This bug resurfaced after implementing in concatenated vector format!

		if epoch != 0 {
			fmt.Println("Epoch ", epoch, "has environments", popstart.NovEnvs)
		}

		pop1 := popstart.Evolve(false, ftraj, json_out, "", epochlength, epoch)
		fmt.Println("End of epoch", epoch)

		if epoch == maxepochs { //Export output population; just before epoch change
			//Update to environment just before epoch change
			pop1.AncEnvs = multicell.CopyCues(pop1.NovEnvs)
			pop1.ToJSON(json_out)
			pop1.Dump_Phenotypes(P_Filename, epochlength)
		}
		dtevol := time.Since(tevol)
		fmt.Println("Time taken to simulate evolution :", dtevol)

		popstart = pop1 //Update population after evolution.
		//fmt.Println("Novel environment before :", popstart.NovEnvs)
		//fmt.Println("Ancestral environment before :", popstart.RefEnvs)

		OldEnvs := multicell.CopyCues(popstart.NovEnvs)
		popstart.AncEnvs = OldEnvs
		popstart.NovEnvs = multicell.ChangeEnvs2(OldEnvs, denv)
		err = multicell.DeepVec3NovTest(popstart.NovEnvs, envtraj)
		if err != nil {
			fmt.Println(err)
		}
		novvec = append(novvec, err == nil)
		//fmt.Println("Novel environment after :", popstart.NovEnvs)
		//fmt.Println("Ancestral environment after :", popstart.AncEnvs)
	}
	err = ftraj.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Trajectory of population written to %s \n", T_Filename)
	fmt.Printf("JSON encoding of evolved population written to %s \n", jfilename)
	fmt.Printf("Phenotypes of population written to %s \n ", P_Filename)
	//fmt.Println("Trajectory of environment :", envtraj)

	fmt.Println("Novelty of environment cue :", novvec)
	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
