package main

import (
	"flag"
	"fmt"
	"log"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string

func main() {
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	gidfilenamePtr := flag.String("geneology", "geneol", "Basename of geneology data files")
	jsoninPtr := flag.String("jsonin", "", "JSON file of input population") //default to empty string
	flag.Parse()

	epochlength := *genPtr
	fmt.Println("epochlength", epochlength)
	Gid_Filename = *gidfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename := fmt.Sprintf("%s_001.json", json_in)
	pop0.ImportPopGz(jfilename)
	multicell.SetParams(pop0.Params)

	fmt.Println("Making DOT genealogy file")
	tdot := time.Now()
	multicell.DOT_Genealogy(Gid_Filename, json_in, epochlength, multicell.GetMaxPop())
	dtdot := time.Since(tdot)
	fmt.Println("Time taken to make dot file :", dtdot)

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
