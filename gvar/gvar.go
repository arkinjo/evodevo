package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var GVar_Filename string //Dump for genotype variance (and entropy?)
//var json_in string       //JSON encoding of initial population; default to empty string
var SS_Anc, SS_Nov float64

func main() {
	log.Println("Starting...")
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	maxdevstepP := flag.Int("maxdevstep", 200, "maximum number of steps for development")
	ngenesP := flag.Int("ngenes", 200, "number of genes")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	//gmodePtr := flag.Int("mode", 1, "Genome to be extracted. 0: Ancestral, 1: Novel") //No need
	//envP := flag.Bool("env", true, "with environment?")
	//ref1Ptr := flag.String("ref1", "", "reference JSON file 1")
	//ref2Ptr := flag.String("ref2", "", "reference JSON file 2")

	//pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	gvarfilenamePtr := flag.String("GVar_file", "gvar", "Filename of genetic variances")
	jsongzinPtr := flag.String("jsongzin", "", "basename of JSON files")

	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	settings.MaxDevStep = *maxdevstepP
	settings.NGenes = *ngenesP
	settings.NCells = *ncellsP

	epochlength := *genPtr

	//refgen1 := *ref1Ptr
	//refgen2 := *ref2Ptr
	//withEnv := *envP

	GVar_Filename = fmt.Sprintf("%s.gvar", *gvarfilenamePtr)

	json_in := *jsongzinPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	/*
		if refgen1 == "" {
			log.Fatal("Must specify JSON reference file 1.")
		}
		if refgen2 == "" {
			log.Fatal("Must specify JSON reference file 2.")
		}
	*/

	//log.Println("epochlength", epochlength)
	fout, err := os.OpenFile(GVar_Filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintf(fout, "#Gen \t Anc \t Nov \n") //Header

	//log.Println("Dumping start \n")
	for gen := 1; gen <= epochlength; gen++ {
		jfilename := fmt.Sprintf("%s_%3.3d.json.gz", json_in, gen)
		pop := multicell.NewPopulation(settings)
		pop.ImportPopGz(jfilename)
		settings = pop.Params
		multicell.SetParams(settings)

		AncPopGVecs := pop.GetFlatGenome(multicell.IAncEnv)
		SSEVec := multicell.GetVarVec(AncPopGVecs)
		SS_Anc = multicell.SumVec(SSEVec)

		NovPopGVecs := pop.GetFlatGenome(multicell.INovEnv)
		SSEVec = multicell.GetVarVec(NovPopGVecs)
		SS_Nov = multicell.SumVec(SSEVec)

		fmt.Fprintf(fout, "%d\t%e\t%e\n", gen, SS_Anc, SS_Nov)

	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	dt := time.Since(t0)

	fmt.Println("Total time taken : ", dt)
}

/*
func avesd(data []float64) (float64, float64) {
	ave := 0.0
	v := 0.0
	fn := float64(len(data))
	for _, d := range data {
		ave += d
	}
	ave /= fn
	for _, d := range data {
		dev := d - ave
		v += dev * dev
	}

	sd := math.Sqrt(v / fn)

	return ave, sd
}
*/
