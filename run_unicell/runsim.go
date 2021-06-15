package main

import (
	"flag"
	//"encoding/json"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/unicell"
)


var T_Filename string = "traj.dat"
var P_Filename string = "phenotypes.dat" //Expressed phenotypes of population
var G_Filename string = "genotypes.dat" //Genome of population
var jsonfile string = "json_test.json"

func main() {
	t0 := time.Now()
	seedPtr := flag.Int("seed", 1, "random seed")
    epochPtr := flag.Int("nepoch", 1, "number of epochs")
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", true, "develop with environmental cue")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 1, "magnitude of environmental changes")
	tfilenamePtr := flag.String("tfilename","traj.txt","name of file of trajectories")
	pfilenamePtr := flag.String("pfilename","phenotypes.txt","name of file of phenotypes")
	gfilenamePtr := flag.String("gfilename","genotypes.txt","name of file of genomes")
    flag.Parse()

	unicell.SetSeed(int64(*seedPtr))
	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = *tfilenamePtr
	P_Filename = *pfilenamePtr
	G_Filename = *gfilenamePtr
	unicell.WithCue = *cuePtr
	unicell.Omega = *omegaPtr


	fout, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file
	if err != nil {
		log.Fatal(err)
	}

	fmt.Fprintln(fout, "Epoch \t Generation \t Fitness \t Cue_Plas \t Obs_Plas \t Utility") //header
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
	pop0 := unicell.NewPopulation(unicell.MaxPop)
	popstart := pop0
	popstart.Env = unicell.RandomEnv(unicell.Nenv,0.5)
	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint) 

	envtraj := make([]unicell.Cue,1) //Trajectory of environment cue

	for epoch := 1; epoch <= maxepochs; epoch++ {
		envtraj = append(envtraj, popstart.Env)

		if epoch != 0 {
			fmt.Println("Epoch ",epoch,"has environment",popstart.Env)
		}
		pop1 := unicell.RecEvolve(T_Filename,epochlength, epoch, &popstart)
		fmt.Println("End of epoch", epoch)
		if epoch == maxepochs {	//JSON encoding is much faster than this.
			fmt.Println("Dumping phenotypes and genotypes")
			t_ext := time.Now()
			popstart.Dump_Phenotypes(P_Filename)
			popstart.Dump_Genotypes(G_Filename)
			dt_ext := time.Since(t_ext)
			fmt.Println("Extraction time : ",dt_ext)
		}
		popstart = pop1  //Update population after evolution.
		OldEnv := popstart.Env.CopyCue()
		popstart.RefEnv = OldEnv
		popstart.Env = popstart.Env.ChangeEnv(denv)
	}

	/*
	jsonpop, err := json.Marshal(pop1) //JSON encoding of population
	if err != nil {
		log.Fatal(err)
	}
	jsonout, err := os.OpenFile(jsonfile, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file
	if err != nil {
		log.Fatal(err)
	}
		
	jsonout.Write(jsonpop)
	if err != nil {
		log.Fatal(err)
		}
	*/

	fmt.Println("Trajectory of population written to",T_Filename)
	fmt.Println("JSON encoding of evolved population written to ",jsonfile)
	fmt.Println("Trajectory of environment :", envtraj)
	
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}
