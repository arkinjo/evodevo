package main

import (
	"encoding/json"
	"flag"
	"fmt"

	//"io/ioutil"
	"log"
	//"math"
	//"os"
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

var JSONequalGs1 float64 //bugtesting variable
var JSONequalGs2 float64 //bugtesting variable

//var test bool = false //false : training mode, true : testing mode
//This file is for playing around; testing behavior of various functions in Go etc.

/*
Objectives:
Fix genome projection bug
Test json population import/export with(out) flag O_Append

*/

func main() {
	//var mu0, mu1 float64
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

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop) //Not cleared
	jsonpop1, err := json.Marshal(pop0)
	if err != nil {
		log.Fatal(err)
	}

	pop1 := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop)
	fmt.Println("Before clearing")
	for _, indiv := range pop1.Indivs {
		fmt.Println(indiv.Genome)
	}
	pop1.ClearGenome()
	fmt.Println("After Clearing")
	for _, indiv := range pop1.Indivs {
		fmt.Println(indiv.Genome)
	}

	err = json.Unmarshal(jsonpop1, &pop1)
	if err != nil {
		log.Fatal(err)
	}

	JSONequalGs1 = multicell.TestEqualGenomes(pop0.Indivs[1].Genome, pop1.Indivs[1].Genome)
	fmt.Println("0 VS 1 :", JSONequalGs1)

	jsonpop2, err := json.Marshal(pop1)
	if err != nil {
		log.Fatal(err)
	}

	pop2 := multicell.NewPopulation(multicell.GetNcells(), multicell.MaxPop)
	pop2.ClearGenome()
	err = json.Unmarshal(jsonpop2, &pop2)
	if err != nil {
		log.Fatal(err)
	}
	JSONequalGs2 = multicell.TestEqualGenomes(pop1.Indivs[1].Genome, pop2.Indivs[1].Genome)
	JSONequalGs = multicell.TestEqualGenomes(pop0.Indivs[1].Genome, pop2.Indivs[1].Genome)
	fmt.Println("1 VS 2 :", JSONequalGs2)
	fmt.Println("0 VS 2 :", JSONequalGs)
	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)

}
