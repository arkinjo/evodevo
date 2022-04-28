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
	state0 := flag.String("state0", "P", "State 0 (one of E, F, G, H, P)")
	state1 := flag.String("state1", "P", "State 1 (one of E, F, G, H, P)")
	ienv0 := flag.Int("ienv0", 0, "0 = Ancestral; 1 = Novel environment")
	ienv1 := flag.Int("ienv1", 1, "0 = Ancestral; 1 = Novel environment")
	flag.Parse()

	multicell.SetMaxPop(*maxpopsizePtr)
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *phenofeedbackPtr, *epigPtr, *HOCPtr)

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

	fmt.Println("Cross-covariance:", *state0, *ienv0, " vs ", *state1, *ienv1)
	mstate0, mstate1, ccmat := pop0.GetCellCrossCov(*state0, *ienv0, *state1, *ienv1)
	fmt.Println("mean", *state0, *ienv0, len(mstate0), ":", mstate0, "\n")
	fmt.Println("mean", *state1, *ienv1, len(mstate1), ":", mstate1, "\n")
	fmt.Println("ccmat[0][0]", ccmat[0][0])
	runSVD(ccmat)
	//	myEigenSym(ccmat)
}

func getTest3by3Matrix() multicell.Dmat {

	dmat := multicell.NewDmat(3, 3)
	dmat[0][0] = 1
	dmat[0][1] = 2
	dmat[0][2] = -3
	dmat[1][0] = 2
	dmat[1][1] = 3
	dmat[1][2] = 4
	dmat[2][0] = -3
	dmat[2][1] = 4
	dmat[2][2] = 5

	return dmat
}

func runSVD(ccmat multicell.Dmat) {
	dim0 := len(ccmat)
	dim1 := len(ccmat[0])
	C := mat.NewDense(dim0, dim1, nil)
	for i, ci := range ccmat {
		for j, v := range ci {
			C.Set(i, j, v)
		}
	}

	var svd mat.SVD
	ok := svd.Factorize(C, mat.SVDFull)
	if !ok {
		log.Fatal("SVD failed.")
	}
	U := mat.NewDense(dim0, dim0, nil)
	V := mat.NewDense(dim1, dim1, nil)
	svd.UTo(U)
	svd.VTo(V)
	vals := svd.Values(nil)
	fmt.Println("\nSValues:", vals)

	for i := 0; i < dim0; i++ {
		fmt.Printf("U\t%f\n", U.At(i, 0))
	}
	for i := 0; i < dim1; i++ {
		fmt.Printf("V\t%f\n", V.At(i, 0))
	}
}

func myEigenSym(m multicell.Dmat) {
	dim := len(m)
	a := mat.NewSymDense(dim, nil)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			a.SetSym(i, j, m[i][j])
		}
	}

	var eig mat.EigenSym
	ok := eig.Factorize(a, true)
	if !ok {
		log.Fatal("EigenSym failed.")
	}
	fmt.Println("EigenValues: ", eig.Values(nil))
	U := mat.NewDense(dim, dim, nil)
	eig.VectorsTo(U)
	for i := 0; i < dim; i++ {
		fmt.Printf("EV\t%f\n", U.At(i, dim-1))
	}
}
