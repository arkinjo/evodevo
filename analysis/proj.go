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
var jfilename string
var json_in string //JSON encoding of initial population; default to empty string
var test bool = false //false : training mode, true : testing mode

func main() {
	t0 := time.Now()
	ncelltypesPtr := flag.Int("celltypes",1,"number of cell types/phenotypes simultaneously trained") //default to unicellular case
    genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	cuePtr := flag.Bool("withCue", false, "develop with environmental cue")
	epigPtr := flag.Bool("epig",false,"Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC",false,"Add layer representing higher order complexes")
	HOIPtr := flag.Bool("HOI",false,"Allow interactions between higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	pgfilenamePtr := flag.String("pgfilename","pg","name of file of projected phenotypes and genotypes") 
	gidfilenamePtr := flag.String("gidfilename","gid","name of file of geneology of ids") 
	jsoninPtr := flag.String("jsonin","","json file of initial population") //default to empty string
    flag.Parse()

	epochlength := *genPtr
	PG_Filename = *pgfilenamePtr
	Gid_Filename = *gidfilenamePtr
	json_in = *jsoninPtr
	multicell.WithCue = *cuePtr
	multicell.Epig = *epigPtr
	multicell.HOC = *HOCPtr
	multicell.HOI = *HOIPtr
	multicell.Omega = *omegaPtr
	multicell.Ncells = *ncelltypesPtr

	pop0 := multicell.NewPopulation(multicell.Ncells,multicell.MaxPop)
	pop1 := multicell.NewPopulation(multicell.Ncells,multicell.MaxPop)
	
	jfilename = fmt.Sprintf("%s_1.json",json_in)
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
	fmt.Println("Successfully imported initial population")
	jfilename = fmt.Sprintf("%s_%d.json",json_in,epochlength)
	fmt.Printf("Importing evolved population from %s \n",jfilename)
	popin, err = os.Open(jfilename)
	if err != nil {
		log.Fatal(err)
	}

	byteValue, _ = ioutil.ReadAll(popin)
	err = json.Unmarshal(byteValue, &pop1)
	if err != nil {
		log.Fatal(err)
	}

	err = popin.Close()
	if err != nil{
		log.Fatal(err)
	}
	fmt.Println("Successfully imported evolved population")

	fmt.Println("Dumping projections")
	tdump := time.Now()
	pop := multicell.NewPopulation(multicell.Ncells,multicell.MaxPop)

	g0 := pop0.GetMeanGenome() //Average genome before evolution
	g1 := pop1.GetMeanGenome() //Average genome after evolution
	Gaxis := multicell.NewGenome()
	multicell.DiffGenomes(Gaxis,g1,g0)
	Gaxis = Gaxis.NormalizeGenome()
	//pop.Envs = pop1.Envs
	//pop.RefEnvs = pop1.RefEnvs

	for gen := 1; gen<=epochlength; gen++ {
		jfilename = fmt.Sprintf("%s_%d.json",json_in,gen)
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
	nanctraj := multicell.DOT_Genealogy(Gid_Filename,json_in,epochlength,multicell.MaxPop)
	fmt.Println(nanctraj)
	dtdot := time.Since(tdot)
	fmt.Println("Time taken to make dot file :",dtdot)
	fmt.Println("Dumping number of ancestors")
	nancfilename = fmt.Sprintf("%s_nanc.dat",Gid_Filename)
	fout, err := os.OpenFile(nancfilename, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644) //create file for recording trajectory
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

	fmt.Printf("Projections written to %s.dat \n",PG_Filename)
	fmt.Printf("Genealogy of final generation written to %s.dot\n",Gid_Filename)
	fmt.Printf("Number of ancestors of final generation written to %s \n",nancfilename)

	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}