package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
	"log"
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

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)

	genome := pop.GetFlatGenome()

	mp0, mg, covp0G := multicell.GetCrossCov(p0, genome, true, true)
	mp1, _, covp1G := multicell.GetCrossCov(p1, genome, true, true)
	dp := multicell.NewVec(lenE)
	multicell.DiffVecs(dp, mp1, mp0)
	lenG := len(mg)

	dirE := mat.NewVecDense(lenE, multicell.CopyVec(denv))
	dirE.ScaleVec(1.0/dirE.Norm(2), dirE)

	dirP := mat.NewVecDense(lenE, multicell.CopyVec(dp))
	dirP.ScaleVec(1.0/dirP.Norm(2), dirP)

	dirG := mat.NewVecDense(lenG, multicell.NewVec(lenG))
	dirG.ScaleVec(1.0/dirG.Norm(2), dirG)

	U0, val0, V0 := multicell.GetSVD(covp0G)
	U1, val1, V1 := multicell.GetSVD(covp1G)

	multicell.ProjectSVD("<Dp0DG>", dirP, dirG, p0, genome, mp0, mg, covp0G, val0, U0, V0)
	multicell.ProjectSVD("<Dp1DG>", dirP, dirG, p1, genome, mp1, mg, covp1G, val1, U1, V1)

	fmt.Println("#V0.V1")
	for i := range val0 {
		u0 := U0.ColView(i)
		v0 := V0.ColView(i)
		for j := range val1 {
			u1 := U1.ColView(j)
			v1 := V1.ColView(j)
			dotu := mat.Dot(u0, u1)
			dotv := mat.Dot(v0, v1)
			fmt.Printf("DotUV\t%d\t%d\t%f\t%f\n", i, j, dotu, dotv)
		}
	}
	fmt.Println("#Uvec")
	for i := 0; i < lenE; i++ {
		fmt.Printf("Uvec01\t%d", i)
		for j := 0; j < 5; j++ {
			fmt.Printf("\t%e\t%e", U0.At(i, j), U1.At(i, j))
		}
		fmt.Println("")
	}
	fmt.Println("#Vvec")
	for i := 0; i < lenG; i++ {
		fmt.Printf("Vvec01\t%d", i)
		for j := 0; j < 5; j++ {
			fmt.Printf("\t%e\t%e", V0.At(i, j), V1.At(i, j))
		}
		fmt.Println("")
	}

}
