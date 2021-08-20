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
var PG_Filename string //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"
//var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	seedPtr := flag.Int("seed", 1, "random seed")
    epochPtr := flag.Int("nepoch", 1, "number of epochs")
	ncelltypesPtr := flag.Int("celltypes",1,"number of cell types/phenotypes simultaneously trained") //default to unicellular case
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", false, "develop with environmental cue")
	epigPtr := flag.Bool("epig",false,"Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC",false,"Add layer representing higher order complexes")
	HOIPtr := flag.Bool("HOI",false,"Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	denvPtr := flag.Int("denv", 2, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename","traj","name of file of trajectories")
	gidfilenamePtr := flag.String("gidfilename","","name of file of geneology of ids") //default to empty string
	jsoninPtr := flag.String("jsonin","","json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout","popout","json file of output population")
	//testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
    flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	maxepochs := *epochPtr
	epochlength := *genPtr
	denv := *denvPtr
	T_Filename = fmt.Sprintf("../analysis/%s.dat",*tfilenamePtr)
	Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.WithCue = *cuePtr
	multicell.Epig = *epigPtr
	multicell.HOC = *HOCPtr
	multicell.HOI = *HOIPtr
	multicell.Omega = *omegaPtr
	multicell.Ncells = *ncelltypesPtr
	//test = *testPtr

	pop0 := multicell.NewPopulation(multicell.Ncells,multicell.MaxPop)
	
	if  json_in != "" { //read input population as a json file, if given
		jfilename := fmt.Sprintf("%s.json",json_in) //Need to be in same folder/directory
		fmt.Printf("Importing initial population from %s \n",jfilename)
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
		gidfilename := fmt.Sprintf("%s_full",Gid_Filename)

		pop1 := multicell.Evolve(true,T_Filename,json_out,gidfilename,epochlength, epoch, &popstart)
		fmt.Println("End of epoch", epoch)

		dtevol := time.Since(tevol)
		fmt.Println("Time taken to simulate evolution :",dtevol)
		popstart.Envs = pop1.Envs.ChangeEnvs(denv) //no updating population needed for test mode
	}

	fmt.Println("Trajectory of population written to",T_Filename)
	fmt.Printf("Genealogy of final generation written to %s.dot\n",Gid_Filename)
	fmt.Printf("JSON encoding of evolved population written to %s.json \n", json_out)
	fmt.Println("Trajectory of environment :", envtraj)
	
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}