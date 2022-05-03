package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"

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

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)

	genome := pop.GetFlatGenome()

	//	mp0, me0, mg, cov0 := multicell.GetCrossCov3(p0, e0, genome, true, true, true)
	//	mp1, me1, _, cov1 := multicell.GetCrossCov3(p1, e1, genome, true, true, true)
	log.Println("GetCrossCov3 .. 0")
	_, _, _, cov0 := multicell.GetCrossCov3(p0, e0, genome, true, true, true)
	log.Println("GetCrossCov3 .. 1")
	/* _, _, _, cov1 := */ multicell.GetCrossCov3(p1, e1, genome, true, true, true)

	s0, _, _, _ := multicell.GetST_HOSVD(cov0)
	for i, s0i := range s0 {
		for j, s0ij := range s0i {
			for k, s0ijk := range s0ij {
				fmt.Println("HOSVD0", i, j, k, s0ijk)
			}
		}
	}

	// dp := multicell.NewVec(lenE)
	// multicell.DiffVecs(dp, mp1, mp0)
	// lenG := len(mg)

	// dirE := mat.NewVecDense(lenE, multicell.CopyVec(denv))
	// dirE.ScaleVec(1.0/dirE.Norm(2), dirE)

	// dirP := mat.NewVecDense(lenE, multicell.CopyVec(dp))
	// dirP.ScaleVec(1.0/dirP.Norm(2), dirP)

	// dirG := mat.NewVecDense(lenG, multicell.NewVec(lenG))
	// dirG.ScaleVec(1.0/dirG.Norm(2), dirG)

}
