package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"

	//	"fmt"
	"github.com/arkinjo/evodevo/multicell"
)

var json_in string
var dp float64

func main() {
	seedPtr := flag.Int("seed", 1, "random seed")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case
	cuestrengthPtr := flag.Float64("cuestrength", 1.0, "control size of var contribution of environmental cue")
	hoistrengthPtr := flag.Float64("hoistrength", 1.0, "control size of var contribution of higher order interactions")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	//omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	//nsamplePtr := flag.Int("nsample", 100, "number of samples")
	flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	//multicell.Omega = *omegaPtr
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *epigPtr, *HOCPtr)
	json_in = *jsoninPtr

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop()) //with randomized genome to start

	if json_in != "" { //read input population as a json file, if given
		fmt.Println("Importing initial population")
		jfilename := fmt.Sprintf("../pops/%s.json", json_in)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}

		byteValue, _ := ioutil.ReadAll(popin)
		pop0.ClearGenome() //Clear genome before unmarshalling json
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

	env0 := multicell.CopyCues(pop0.RefEnvs)

	lambda2 := 0.0

	for _, indiv := range pop0.Indivs {
		//indiv0 := multicell.NewIndiv(i)
		conv0 := false
		conv1 := false
		//indiv0.Genome.Randomize()
		_, err0 := indiv.Copies[multicell.IAncEnv].DevCells(indiv.Genome, env0)
		_, err1 := indiv.Copies[multicell.INovEnv].DevCells(indiv.Genome, env0)
		if err0 == nil {
			conv0 = true
		}
		if err1 == nil {
			conv1 = true
		}
		for j, cell := range indiv.Copies[multicell.IAncEnv].Ctypes { //Loop over all cells
			dp = multicell.Dist2Vecs(cell.P, indiv.Copies[multicell.INovEnv].Ctypes[j].P)
			if dp > 1.0e-2 {
				fmt.Println("ID:", indiv.Id, "Diff:", dp, conv0, conv1)
			}
		}
		lambda2 += dp
	}
	fmt.Println("Population Instability Measure:", 0.5*math.Log(lambda2))
}
