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
var jsongz_in string //gzip compressed JSON encoding of initial population; default to empty string
var jsongz_out string = "popout"
var jfilename string

func main() {
	t0 := time.Now()
	testP := flag.Bool("test", false, "Test run or not")
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	maxdevstepP := flag.Int("maxdevstep", 200, "maximum number of steps for development")
	ngenesP := flag.Int("ngenes", 200, "Number of genes")
	nenvP := flag.Int("nenv", 200, "Number of environmental cues/traits")
	nselP := flag.Int("nsel", 40, "Number of environmental cues/traits for selection")
	ncellsP := flag.Int("ncells", 1, "Number of cell types")
	withcueP := flag.Bool("cue", true, "With environmental cue")
	flayerP := flag.Bool("layerF", true, "Epigenetic layer")
	hlayerP := flag.Bool("layerH", true, "Higher order complexes")
	jlayerP := flag.Bool("layerJ", true, "Interactions in higher order interactions")
	pfbackP := flag.Bool("pfback", true, "Phenotype feedback to input")

	noiseP := flag.Float64("noise", 0.05, "Strength of environmental noise")
	mutP := flag.Float64("mut", 0.005, "Mutation rate")
	tauFP := flag.Float64("tauF", 0.2, "Decay rate of the f layer")
	denEP := flag.Float64("dE", 0.02, "Density of E")
	denFP := flag.Float64("dF", 0.02, "Density of F")
	denGP := flag.Float64("dG", 0.02, "Density of G")
	denHP := flag.Float64("dH", 0.02, "Density of H")
	denJP := flag.Float64("dJ", 0.02, "Density of J")
	denPP := flag.Float64("dP", 0.02, "Density of P")

	seedPtr := flag.Int("seed", 13, "random seed")
	seed_cuePtr := flag.Int("seed_cue", 7, "random seed for environmental cue")
	epochPtr := flag.Int("nepoch", 20, "number of epochs")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")

	denvPtr := flag.Int("denv", 100, "magnitude of environmental change")
	tfilenamePtr := flag.String("traj_file", "traj.dat", "filename of trajectories")
	jsongzinPtr := flag.String("jsongzin", "", "json file of input population") //default to empty string
	jsongzoutPtr := flag.String("jsongzout", "popout", "json file of output population")
	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	settings.MaxDevStep = *maxdevstepP
	settings.NGenes = *ngenesP
	settings.NEnv = *nenvP
	settings.NSel = *nselP
	settings.NCells = *ncellsP
	settings.WithCue = *withcueP
	settings.WithCue = *withcueP
	settings.FLayer = *flayerP
	settings.HLayer = *hlayerP
	settings.JLayer = *jlayerP
	settings.Pfback = *pfbackP
	settings.SDNoise = *noiseP
	settings.MutRate = *mutP
	settings.TauF = *tauFP
	settings.DensityE = *denEP
	settings.DensityF = *denFP
	settings.DensityG = *denGP
	settings.DensityH = *denHP
	settings.DensityJ = *denJP
	settings.DensityP = *denPP

	log.Println("seed=", *seedPtr, "seed_cue=", *seed_cuePtr)
	multicell.SetSeed(int64(*seedPtr))
	multicell.SetSeedCue(int64(*seed_cuePtr))

	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = *tfilenamePtr
	jsongz_in = *jsongzinPtr
	jsongz_out = *jsongzoutPtr
	test_flag := *testP

	pop0 := multicell.NewPopulation(settings)

	if jsongz_in != "" { //read input population as a .json.gz file, if given
		pop0.ImportPopGz(jsongz_in)
	}

	pop0.Params.SDNoise = settings.SDNoise
	pop0.Params.MutRate = settings.MutRate
	multicell.SetParams(pop0.Params)
	if jsongz_in == "" {
		fmt.Println("Randomizing initial population")
		pop0.RandomizeGenome()
	}

	ftraj, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0
	if jsongz_in != "" {
		popstart.ChangeEnvs(denv)
	} else {
		popstart.SetRandomNovEnvs()
	}

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

		pop1 := popstart.Evolve(test_flag, ftraj, jsongz_out, epochlength, epoch)
		fmt.Println("End of epoch", epoch)

		if !test_flag && epoch == maxepochs { //Export output population; just before epoch change
			pop1.ExportPopGz(jsongz_out)
		}
		dtevol := time.Since(tevol)
		fmt.Println("Time taken to simulate evolution :", dtevol)

		popstart = pop1 //Update population after evolution.
		popstart.ChangeEnvs(denv)
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

