package main

// practice using gonum

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
	"io/ioutil"
	"log"
	//	"math/rand"
	"os"
)

func main() {
	maxpopsizePtr := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular c
	cuestrengthPtr := flag.Float64("cuestrength", 1.0, "control size of var contribution of environmental cue")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	phenofeedbackPtr := flag.Bool("pheno_feedback", false, "controls phenotype feedback into regulation")
	hoistrengthPtr := flag.Float64("hoistrength", 1.0, "control size of var contribution of higher order interactions")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	jsoninPtr := flag.String("jsonin", "", "json file of input population") //default to empty string
	flag.Parse()

	multicell.SetMaxPop(*maxpopsizePtr)
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *phenofeedbackPtr, *epigPtr, *HOCPtr)

	ncells := multicell.GetNcells()
	nenv := multicell.GetNenv()

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop())
	pop0.ClearGenome()

	fmt.Println("Importing population from ", *jsoninPtr)
	popin, err := os.Open(*jsoninPtr)
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

	ccdim := ncells + nenv
	mphen, menv, ccmat := pop0.GetPhenoEnvCC(multicell.IAncEnv)
	fmt.Println("mean phenotype", mphen)
	fmt.Println("mean envs", menv)
	C := mat.NewDense(ccdim, ccdim, nil)
	for i := range ccmat {
		for j, v := range ccmat[0] {
			C.Set(i, j, v)
		}
	}
	var svd mat.SVD
	ok := svd.Factorize(C, mat.SVDFull)
	if !ok {
		log.Fatal("SVD failed.")
	}
	U := mat.NewDense(ccdim, ccdim, nil)
	V := mat.NewDense(ccdim, ccdim, nil)
	svd.UTo(U)
	svd.VTo(V)
	vals := svd.Values(nil)
	fmt.Println("Values:", vals)

	for i := 0; i < ccdim; i++ {
		fmt.Printf("SVec\t%f\t%f:", U.At(i, 0), V.At(i, 0))
		fmt.Println()

	}
}
