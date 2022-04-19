package main

import (
	"encoding/json"
	//"errors"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var PG_Filename string  //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string
var P_Filename string

//var json_out string

var CopyequalGs float64 //bugtesting variable
var JSONequalGs float64 //bugtesting variable
var DevequalGs float64  //bugtesting variable

//var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	//seedPtr := flag.Int("seed", 1, "random seed")
	//epochPtr := flag.Int("nepoch", 1, "number of epochs")
	maxpopsizePtr := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case
	//genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	refgenPtr := flag.Int("refgen", 50, "reference generation for evolved genotype")
	cuestrengthPtr := flag.Float64("cuestrength", 1.0, "control size of var contribution of environmental cue")
	phenofeedbackPtr := flag.Bool("pheno_feedback", false, "controls phenotype feedback into regulation")
	hoistrengthPtr := flag.Float64("hoistrength", 1.0, "control size of var contribution of higher order interactions")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	//HOIPtr := flag.Bool("HOI", true, "Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	//denvPtr := flag.Int("denv", 10, "magnitude of environmental change")
	//tfilenamePtr := flag.String("tfilename", "traj", "name of file of trajectories")
	//pgfilenamePtr := flag.String("pgfilename", "", "name of file of projected phenotypes and genotypes") //default to empty string
	//gidfilenamePtr := flag.String("gidfilename", "", "name of file of geneology of ids")                 //default to empty string
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	pfilenamePtr := flag.String("pfilename", "pvec", "name of file of phenotypes")
	//jsonoutPtr := flag.String("jsonout", "popout", "json file of output population")
	//testPtr := flag.Bool("test",false,"test mode if true, defaults to train mode")
	flag.Parse()

	//multicell.SetSeed(int64(*seedPtr))
	//maxepochs := *epochPtr
	//epochlength := *genPtr
	refgen := *refgenPtr
	//denv := *denvPtr
	//T_Filename = fmt.Sprintf("../analysis/%s.dat", *tfilenamePtr)
	//PG_Filename = *pgfilenamePtr
	//Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	P_Filename = fmt.Sprintf("../analysis/%s_%d.dat", *pfilenamePtr, refgen) //all phenotypes go to analysis directory

	//json_out = *jsonoutPtr
	multicell.Omega = *omegaPtr

	multicell.SetMaxPop(*maxpopsizePtr)
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *phenofeedbackPtr, *epigPtr, *HOCPtr)

	pop := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop()) //with randomized genome to start
	//pop1 := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop()) //with randomized genome to start

	fmt.Println("Importing population")
	jfilename := fmt.Sprintf("../pops/%s.json", json_in)
	popin, err := os.Open(jfilename)
	if err != nil {
		log.Fatal(err)
	}

	byteValue, _ := ioutil.ReadAll(popin)
	pop.ClearGenome() //Clear genome before unmarshalling json
	err = json.Unmarshal(byteValue, &pop)
	if err != nil {
		log.Fatal(err)
	}

	err = popin.Close()
	if err != nil {
		log.Fatal(err)
	}

	pop.Dump_Phenotypes(P_Filename, refgen)

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
