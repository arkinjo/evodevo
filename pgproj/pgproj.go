package main

import (
	"flag"
	"fmt"
	"log"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var PG_Filename string //Dump for phenotypes and genotypes
var json_in string     //JSON encoding of initial population; default to empty string

func main() {
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case

	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	refgenPtr := flag.Int("refgen", 50, "reference generation for evolved genotype")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "JSON file of input population") //default to empty string

	flag.Parse()

	epochlength := *genPtr
	fmt.Println("epochlength", epochlength)
	refgen := *refgenPtr
	PG_Filename = *pgfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename := fmt.Sprintf("%s_001.json", json_in)
	pop0.FromJSON(jfilename)
	multicell.SetParams(pop0.Params)

	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename = fmt.Sprintf("%s_%3.3d.json", json_in, refgen)
	fmt.Println("Reference population :", jfilename)
	pop1.FromJSON(jfilename)

	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	fmt.Println("Dumping projections")
	tdump := time.Now()

	g0 := pop0.GetMeanGenome()
	g1 := pop1.GetMeanGenome()
	Gaxis := multicell.NewGenome()
	multicell.DiffGenomes(&Gaxis, &g1, &g0) //Pointers for bugfix; don't ask why; it just works!
	Gaxis = Gaxis.NormalizeGenome()
	Paxis := pop1.Get_Environment_Axis() //Measure everything in direction of ancestral -> novel environment

	for gen := 1; gen <= epochlength; gen++ { //Also project population after pulling back to ancestral environment.
		jfilename := fmt.Sprintf("%s_%3.3d.json", json_in, gen)
		pop := multicell.NewPopulation(*ncellsP, *maxpopP)
		pop.FromJSON(jfilename)
		pop.Dump_Projections(PG_Filename, gen, Gaxis, Paxis)
	}

	dtdump := time.Since(tdump)
	fmt.Println("Time taken to dump projections :", dtdump)
	fmt.Printf("Projections written to %s_*.dat \n", PG_Filename)
	dt := time.Since(t0)

	fmt.Println("Total time taken : ", dt)
}
