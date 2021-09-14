package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var PG_Filename string  //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"

var CopyequalGs float64 //bugtesting variable
var JSONequalGs float64 //bugtesting variable
var DevequalGs float64  //bugtesting variable
var Gdiff float64

//var test bool = false //false : training mode, true : testing mode
//This file is for playing around; testing behavior of various functions in Go etc.

/*
Objectives:
Fix genome projection bug
Test json population import/export with(out) flag O_Append

*/

func main() {
	var mu0, mu1 float64
	t0 := time.Now()
	fmt.Println("Hello, world!")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case
	cuePtr := flag.Bool("withCue", true, "develop with environmental cue")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	HOIPtr := flag.Bool("HOI", true, "Allow interactions between higher order complexes")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")

	json_in = *jsoninPtr
	json_out = *jsonoutPtr

	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuePtr, *epigPtr, *HOCPtr, *HOIPtr)

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop)
	pop1 := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop)

	if json_in != "" { //read input population as a json file, if given
		fmt.Println("Importing initial population")
		jfilename := fmt.Sprintf("../pops/%s.json", json_in)
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
		if err != nil {
			log.Fatal(err)
		}
		fmt.Println("Successfully imported population")
	}

	playenv := multicell.RandomEnvs(multicell.GetNcells(), multicell.GetNenv(), 0.5) //Random environment cue
	pop0.Envs = playenv
	zeroGenome := multicell.NewGenome()

	for _, indiv := range pop0.Indivs {
		Gdiff = multicell.TestEqualGenomes(zeroGenome, indiv.Genome)
		fmt.Println("ID : ", indiv.Id, "; Genome metric value :", Gdiff) //Check whether this is non-zero
		mu0 += Gdiff / float64(1000)
	}
	fmt.Println("Marshalling and unmarshalling")

	jsonpop, err := json.Marshal(pop0) //Marshal json encoding of pop0
	if err != nil {
		log.Fatal(err)
	}
	err = json.Unmarshal(jsonpop, &pop1) //Unmarshal into pop1
	if err != nil {
		log.Fatal(err)
	}

	for _, indiv := range pop1.Indivs {
		Gdiff = multicell.TestEqualGenomes(zeroGenome, indiv.Genome)
		fmt.Println("ID : ", indiv.Id, "; Genome metric value :", Gdiff) //Check whether this is non-zero
		mu1 += Gdiff / float64(1000)
	}

	JSONequalGs = multicell.TestEqualPopGenomes(pop0, pop1)
	fmt.Println("Total Marshaller error :", JSONequalGs)
	fmt.Println("Mean marshaller error :", JSONequalGs/float64(1000))
	fmt.Println("mu_0-mu_1=", mu0-mu1)
	fmt.Println("Error of errors :", JSONequalGs/float64(1000)-math.Abs(mu0-mu1))

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)

}
