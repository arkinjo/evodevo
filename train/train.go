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
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"

func main() {
	t0 := time.Now()
	seedPtr := flag.Int("seed", 1, "random seed")
    epochPtr := flag.Int("nepoch", 20, "number of epochs")
	ncelltypesPtr := flag.Int("celltypes",1,"number of cell types/phenotypes simultaneously trained") //default to unicellular case
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", false, "develop with environmental cue")
	epigPtr := flag.Bool("epig",false,"Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC",false,"Add layer representing higher order complexes")
	HOIPtr := flag.Bool("HOI",false,"Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 20, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename","traj","name of file of trajectories")
	jsoninPtr := flag.String("jsonin","","json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout","popout","json file of output population")
    flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = fmt.Sprintf("%s.dat",*tfilenamePtr)
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.WithCue = *cuePtr
	multicell.Epig = *epigPtr
	multicell.HOC = *HOCPtr
	multicell.HOI = *HOIPtr
	multicell.Omega = *omegaPtr
	multicell.Ncells = *ncelltypesPtr

	pop0 := multicell.NewPopulation(multicell.Ncells,multicell.MaxPop)
	
	if  json_in != "" { //read input population as a json file, if given
		jfilename := fmt.Sprintf("%s.json",json_in)
		fmt.Printf("Importing initial population from %s.json \n",jfilename)
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

	fmt.Fprintln(fout, "Epoch \t Generation \t Fitness \t Cue_Plas \t Obs_Plas \t Polyphenism \t Utility") //header
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0
	popstart.Envs = multicell.RandomEnvs(multicell.Ncells,multicell.Nenv,0.5)
	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)
	
	envtraj := make([]multicell.Cues,1) //Trajectory of environment cue
	envtraj[0] = popstart.RefEnvs

	for epoch := 1; epoch <= maxepochs; epoch++ {
		tevol := time.Now()
		envtraj = append(envtraj, popstart.Envs)

		if epoch != 0 {
			fmt.Println("Epoch ",epoch,"has environments",popstart.Envs)
		}

		pop1 := multicell.Evolve(false,T_Filename,json_out,"",epochlength, epoch, &popstart)
		fmt.Println("End of epoch", epoch)

		if epoch == maxepochs { //Export output population
			jfilename := fmt.Sprintf("../test/%s.json",json_out) //export output population to test file
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
		dtevol := time.Since(tevol)
		fmt.Println("Time taken to simulate evolution :",dtevol)

		popstart = pop1  //Update population after evolution.
		OldEnv := popstart.Envs.CopyCues()
		popstart.RefEnvs = OldEnv
		popstart.Envs = OldEnv.ChangeEnv(denv)
	}

	fmt.Println("Trajectory of population written to",T_Filename)
	fmt.Printf("JSON encoding of evolved population written to %s.json \n", json_out)
	fmt.Println("Trajectory of environment :", envtraj)
	
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}