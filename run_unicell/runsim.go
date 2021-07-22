package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/unicell"
)


var T_Filename string = "traj"
var PG_Filename string //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"
var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	seedPtr := flag.Int("seed", 1, "random seed")
    epochPtr := flag.Int("nepoch", 1, "number of epochs")
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", true, "develop with environmental cue")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 2, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename","traj","name of file of trajectories")
	pgfilenamePtr := flag.String("pgfilename","","name of file of projected phenotypes and genotypes") //default to empty string
	gidfilenamePtr := flag.String("gidfilename","","name of file of geneology of ids") //default to empty string
	jsoninPtr := flag.String("jsonin","","json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout","popout","json file of output population")
	testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
    flag.Parse()

	unicell.SetSeed(int64(*seedPtr))
	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = fmt.Sprintf("%s.dat",*tfilenamePtr)
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	unicell.WithCue = *cuePtr
	unicell.Omega = *omegaPtr
	test = *testPtr

	pop0 := unicell.NewPopulation(unicell.MaxPop)
	
	if  json_in != "" { //read input population as a json file, if given
		fmt.Println("Importing initial population")
		jfilename := fmt.Sprintf("%s.json",json_in)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}

		byteValue, _ := ioutil.ReadAll(popin)
		err = json.Unmarshal(byteValue, &pop0)
		if err != nil {
			log.Fatal(err)
		}

		err = popin.Close()
		if err != nil{
			log.Fatal(err)
		}
		fmt.Println("Successfully imported population")
	}

	fout, err := os.OpenFile(T_Filename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}

	fmt.Fprintln(fout, "Epoch \t Generation \t Fitness \t Cue_Plas \t Obs_Plas \t Utility") //header
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0
	popstart.Env = unicell.RandomEnv(unicell.Nenv,0.5)
	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)
	
	envtraj := make([]unicell.Cue,1) //Trajectory of environment cue
	envtraj[0] = popstart.RefEnv

	for epoch := 1; epoch <= maxepochs; epoch++ {
		envtraj = append(envtraj, popstart.Env)

		if epoch != 0 {
			fmt.Println("Epoch ",epoch,"has environment",popstart.Env)
		}

		pop1 := unicell.Evolve(test,T_Filename,json_out,Gid_Filename,epochlength, epoch, &popstart)
		fmt.Println("End of epoch", epoch)

		if epoch == maxepochs { //Export output population
			jfilename := fmt.Sprintf("%s.json",json_out)
			jsonpop, err := json.Marshal(pop1) //JSON encoding of population as byte array
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
		}
		if test { //Randomly generate new environment each epoch in test mode
			popstart.Env = unicell.RandomEnv(unicell.Nenv,0.5)
		} else { //Update population in training mode
			popstart = pop1  //Update population after evolution.
			OldEnv := popstart.Env.CopyCue()
			popstart.RefEnv = OldEnv
			popstart.Env = OldEnv.ChangeEnv(denv)
		}
	}

	fmt.Println("Trajectory of population written to",T_Filename,".dat")
	fmt.Println("Projections written to",PG_Filename,".dat")
	fmt.Println("JSON encoding of evolved population written to ",json_out)
	fmt.Println("Trajectory of environment :", envtraj)
	
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}