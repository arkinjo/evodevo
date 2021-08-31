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
    //epochPtr := flag.Int("nepoch", 1, "number of epochs")
	ncelltypesPtr := flag.Int("celltypes",1,"number of cell types/phenotypes simultaneously trained") //default to unicellular case
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", false, "develop with environmental cue")
	epigPtr := flag.Bool("epig",false,"Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC",false,"Add layer representing higher order complexes")
	HOIPtr := flag.Bool("HOI",false,"Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	//denvPtr := flag.Int("denv", 2, "magnitude of environmental change")
	tfilenamePtr := flag.String("tfilename","traj","name of file of trajectories")
	pgfilenamePtr := flag.String("pgfilename","","name of file of projected phenotypes and genotypes") //default to empty string
	gidfilenamePtr := flag.String("gidfilename","","name of file of geneology of ids") //default to empty string
	jsoninPtr := flag.String("jsonin","","json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout","popout","json file of output population")
	//testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
    flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	//maxepochs := *epochPtr
	epochlength := *genPtr
	//denv := *denvPtr
	T_Filename = fmt.Sprintf("%s.dat",*tfilenamePtr)
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	json_out = *jsonoutPtr
	multicell.WithCue = *cuePtr
	multicell.Epig = *epigPtr
	multicell.HOC = *HOCPtr
	multicell.HOI = *HOIPtr
	multicell.Omega = *omegaPtr
	
	multicell.SetNcells(*ncelltypesPtr)
	//multicell.Ncells = *ncelltypesPtr
	//multicell.SelStrength = 0.25/float64(multicell.Ncells)
	//test = *testPtr

	pop0 := multicell.NewPopulation(multicell.GetNcells(),multicell.MaxPop)
	
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

	fmt.Fprintln(fout, "Epoch \t Generation \t Fitness \t Cue_Plas \t Obs_Plas \t Polyphenism \t Diversity \t Utility") //header
	pop0.Envs = pop0.RefEnvs //Generation zero; just before environment change
	pop0.DevPop(0)

	fmt.Fprintf(fout, "1 \t 0 \t %e \t %e \t %e \t %e \t %e \t %e \n",pop0.GetMeanFitness(),pop0.GetMeanCuePlasticity(),pop0.GetMeanObsPlasticity(),pop0.GetMeanPp(),pop0.GetDiversity(),pop0.GetMeanUtility())

	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
	
	jfilename := fmt.Sprintf("../analysis/%s_0.json",json_out)
	jsonpop, err := json.Marshal(pop0) //JSON encoding of population as byte array
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

	err = popout.Close()
	if err != nil {
		log.Fatal(err)
	}

	popstart := pop0
	popstart.Envs = multicell.RandomEnvs(multicell.GetNcells(),multicell.Nenv,0.5)
	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)
	
	envtraj := make([]multicell.Cues,1) //Trajectory of environment cue
	envtraj[0] = popstart.RefEnvs

	//for epoch := 1; epoch <= maxepochs; epoch++ {
	tevol := time.Now()
	envtraj = append(envtraj, popstart.Envs)

	fmt.Println("Novel environment :",popstart.Envs)
	gidfilename := fmt.Sprintf("%s_full",Gid_Filename)

	pop1 := multicell.Evolve(true,T_Filename,json_out,gidfilename,epochlength, 1, &popstart)
	fmt.Println("End of epoch")

	/*
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
	*/
	dtevol := time.Since(tevol)
	fmt.Println("Time taken to simulate evolution :",dtevol)
	fmt.Println("Dumping projections")
	tdump := time.Now()
	pop := multicell.NewPopulation(multicell.GetNcells(),multicell.MaxPop)
	g0 := pop0.GetMeanGenome()
	g1 := pop1.GetMeanGenome()
	Gaxis := multicell.NewGenome()
	multicell.DiffGenomes(Gaxis,g1,g0)
	Gaxis = Gaxis.NormalizeGenome()

	for gen := 0; gen<=epochlength; gen++ {
		jfilename := fmt.Sprintf("../analysis/%s_%d.json",json_out,gen)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}
				
		byteValue, _ := ioutil.ReadAll(popin)
		err = json.Unmarshal(byteValue, &pop)
		if err != nil {
			log.Fatal(err)
		}

		err = popin.Close()
		if err != nil{
			log.Fatal(err)
		}

		pop.Envs = pop1.Envs
		pop.RefEnvs = pop1.RefEnvs
		pop.Dump_Projections(PG_Filename,gen,Gaxis)
	}
	dtdump := time.Since(tdump)
	fmt.Println("Time taken to dump projections :",dtdump)
	fmt.Println("Making DOT genealogy file")
	tdot := time.Now()
	nanctraj := multicell.DOT_Genealogy(Gid_Filename,json_out,epochlength,multicell.MaxPop)
	fmt.Println(nanctraj)
	dtdot := time.Since(tdot)
	fmt.Println("Time taken to make dot file :",dtdot)
	fmt.Println("Dumping number of ancestors")
	nancfilename = fmt.Sprintf("../analysis/%s_nanc.dat",Gid_Filename)
	fout, err = os.OpenFile(nancfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "Generation \t Ancestors")
	for i, n := range(nanctraj){
		fmt.Fprintf(fout,"%d\t%d\n",i+1,n)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Trajectory of population written to",T_Filename)
	fmt.Printf("Projections written to %s.dat \n",PG_Filename)
	fmt.Printf("Genealogy of final generation written to %s.dot\n",Gid_Filename)
	fmt.Printf("Number of ancestors of final generation written to %s.dat\n",nancfilename)
	fmt.Printf("JSON encoding of evolved population written to %s.json \n", json_out)
	fmt.Println("Trajectory of environment :", envtraj)
	
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}