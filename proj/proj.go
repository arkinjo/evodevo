package main

import (
	"flag"
	"fmt"
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
var json_out string = "popout"

var CopyequalGs float64 //bugtesting variable
var JSONequalGs float64 //bugtesting variable
var DevequalGs float64  //bugtesting variable

//var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case

	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	refgenPtr := flag.Int("refgen", 50, "reference generation for evolved genotype")

	tfilenamePtr := flag.String("tfilename", "traj", "Filename of trajectories")
	pgfilenamePtr := flag.String("pgfilename", "phenogeno", "Filename of projected phenotypes and genotypes") //default to empty string
	gidfilenamePtr := flag.String("gidfilename", "geneol", "Filename of geneology of ids")                    //default to empty string
	jsoninPtr := flag.String("jsonin", "", "JSON file of input population")                                   //default to empty string

	flag.Parse()

	epochlength := *genPtr
	fmt.Println("epochlength", epochlength)
	refgen := *refgenPtr
	T_Filename = fmt.Sprintf("../analysis/%s.dat", *tfilenamePtr)
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename := fmt.Sprintf("../pops/%s_001.json", json_in)
	pop0.FromJSON(jfilename)
	multicell.SetParams(pop0.Params)

	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename = fmt.Sprintf("../pops/%s_%3.3d.json", json_in, refgen)
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
		jfilename := fmt.Sprintf("../pops/%s_%3.3d.json", json_in, gen)
		pop := multicell.NewPopulation(*ncellsP, *maxpopP)
		pop.FromJSON(jfilename)
		pop.Dump_Projections(PG_Filename, gen, Gaxis, Paxis)
	}
	dtdump := time.Since(tdump)
	fmt.Println("Time taken to dump projections :", dtdump)
	fmt.Println("Making DOT genealogy file")
	tdot := time.Now()
	nanctraj := multicell.DOT_Genealogy(Gid_Filename, json_in, epochlength, multicell.GetMaxPop())
	dtdot := time.Since(tdot)
	fmt.Println("Time taken to make dot file :", dtdot)
	fmt.Println("Dumping number of ancestors")
	nancfilename = fmt.Sprintf("../analysis/%s_nanc.dat", Gid_Filename)
	fout, err := os.OpenFile(nancfilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644) //create file for recording trajectory
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "Generation \t Ancestors")
	for i, n := range nanctraj {
		fmt.Fprintf(fout, "%d\t%d\n", i+1, n)
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Projections written to %s_*.dat \n", PG_Filename)
	fmt.Printf("Genealogy of final generation written to %s.dot\n", Gid_Filename)
	fmt.Printf("Number of ancestors of final generation written to %s\n", nancfilename)
	fmt.Printf("JSON encoding of populations written to %s_*.json \n", json_in)

	dt := time.Since(t0)
	fmt.Println("Total time taken : ", dt)
}
