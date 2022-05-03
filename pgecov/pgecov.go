package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"

	"log"
	//	"math"

	"github.com/arkinjo/evodevo/multicell"
	//	"gonum.org/v1/gonum/mat"
)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")

	jsonP := flag.String("jsonin", "", "json file of population")

	flag.Parse()

	pop := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsonP != "" {
		pop.FromJSON(*jsonP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	lenE := len(env0)
	denv := multicell.NewVec(lenE)
	multicell.DiffVecs(denv, env1, env0)

	e0 := pop.GetFlatStateVec("E", 0)
	e1 := pop.GetFlatStateVec("E", 1)
	dele := make([][]float64, 0)
	for i, e := range e0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, e1[i], e)
		dele = append(dele, d)
	}

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)
	delp := make([][]float64, 0)
	for i, p := range p0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, p1[i], p)
		delp = append(delp, d)
	}

	genome := pop.GetFlatGenome()

	_, _, _, cov := multicell.GetCrossCov3(delp, dele, genome, true, true, true)
	svals, _, _, _ := multicell.GetFastHOSVD(cov)
	for _, e := range svals {
		fmt.Printf("HOSVD1 %d %d %d %e\n", e.I, e.J, e.K, e.V)
	}
}
