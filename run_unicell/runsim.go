package main

import (
	"flag"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/unicell"
)

//var epoch int
//var train_jsonfile string = "train.txt"
//var test_jsonfile string = "test.txt"
var jsonfile string = "test.json"

func main() {
	t0 := time.Now()
        seedPtr := flag.Int("seed", 11, "random seed")
        epochPtr := flag.Int("nepoch", 1, "number of epochs")
        genPtr := flag.Int("ngen", 100, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", true, "develop with environmental cue")
	//        denvPtr := flag.Int("denv", 1, "magnitude of environmental changes")
        flag.Parse()

	unicell.SetSeed(int64(*seedPtr))
	maxepochs := *epochPtr
	epochlength := *genPtr
	//	denv := *denvPtr

	unicell.WithCue = *cuePtr

	fout, err := os.OpenFile(unicell.Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "Epoch", "Generation", "Fitness", "Plasticity", "Utility")
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
	pop0 := unicell.NewPopulation(unicell.MaxPop)
	popstart := pop0
	popstart.Env = unicell.RandomBoolEnv(unicell.Nenv,0.5)

	envtraj := make([]unicell.Cue,1) //Trajectory of environment cue

	for epoch := 1; epoch <= maxepochs; epoch++ {
		envtraj = append(envtraj, popstart.Env)

		if epoch != 0 {
			fmt.Println("Epoch ",epoch,"has environment",popstart.Env)
		}
		
		pop1 := unicell.RecEvolve(epochlength, &popstart, epoch)
		fmt.Println("End of epoch", epoch)
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
		//fmt.Println("JSON encoding of evolved population written to ",jsonfile)
		
		popstart = pop1  //Update population after evolution.
		OldEnv := popstart.Env.CopyCue()
		popstart.RefEnv = OldEnv
		//popstart.Env = unicell.RandomEnv(unicell.Nenv, 0.5) //Change environment at end of epoch.
		//fmt.Println("Trajectory of environment :", envtraj)

		popstart.Env = OldEnv.ChangeBoolEnv(1) 
	}

	fmt.Println("Trajectory of population written to",unicell.Filename)
	fmt.Println("JSON encoding of evolved population written to ",jsonfile)
	fmt.Println("Trajectory of environment :", envtraj)
	
	/*
	jsonin := make([]byte,n)
	jsonout.Read(jsonin)
	json.Unmarshal(jsonin,&pop0)
	fmt.Println("Before :",pop1)
	fmt.Println("After :",pop0)
	*/

	dt := time.Since(t0)
	fmt.Println("Total time taken = ",dt)
}
