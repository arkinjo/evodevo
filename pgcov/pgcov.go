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

	json0P := flag.String("json0", "", "json file of initial population")
	json1P := flag.String("json1", "", "json file of final population")

	flag.Parse()

	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *json0P != "" {
		pop0.FromJSON(*json0P)
		multicell.SetParams(pop0.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -json0=filename.")
	}
	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *json1P != "" {
		pop1.FromJSON(*json1P)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the reference JSON file with -json1=filename.")
	}

	env0 := multicell.FlattenEnvs(pop0.AncEnvs)
	env1 := multicell.FlattenEnvs(pop1.NovEnvs)
	dim := len(env0)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)

	e0 := pop0.GetFlatStateVec("E", 0)
	e1 := pop1.GetFlatStateVec("E", 1)
	me0 := multicell.GetMeanVec(e0)
	me1 := multicell.GetMeanVec(e1)
	de := multicell.NewVec(dim)
	multicell.DiffVecs(de, me1, me0)

	p0 := pop0.GetFlatStateVec("P", 0)
	p1 := pop1.GetFlatStateVec("P", 1)
	mp0 := multicell.GetMeanVec(p0)
	mp1 := multicell.GetMeanVec(p1)
	dp := multicell.NewVec(dim)
	multicell.DiffVecs(dp, mp1, mp0)

	G0 := pop0.GetFlatGenome()
	G1 := pop1.GetFlatGenome()
	mg0 := multicell.GetMeanVec(G0)
	mg1 := multicell.GetMeanVec(G1)
	lenG := len(mg0)
	dg := multicell.NewVec(lenG)
	multicell.DiffVecs(dg, mg1, mg0)

	deltaP := make([]multicell.Vec, 0)
	deltaE := make([]multicell.Vec, 0)
	deltaG := make([]multicell.Vec, 0)
	for k, p := range p1 {
		d := multicell.NewVec(dim)
		multicell.DiffVecs(d, p, p0[k])
		deltaP = append(deltaP, d)

		e := multicell.NewVec(dim)
		multicell.DiffVecs(e, e1[k], e0[k])
		deltaE = append(deltaE, e)

		G := multicell.NewVec(lenG)
		multicell.DiffVecs(G, G1[k], G0[k])
		deltaG = append(deltaG, G)
	}

	dirE := mat.NewVecDense(dim, multicell.CopyVec(denv))
	dirE.ScaleVec(1.0/dirE.Norm(2), dirE)

	dirP := mat.NewVecDense(dim, multicell.CopyVec(dp))
	dirP.ScaleVec(1.0/dirP.Norm(2), dirP)

	dirG := mat.NewVecDense(len(mg0), dg)
	dirG.ScaleVec(1.0/dirG.Norm(2), dirG)

	pgfn2 := 0.0
	for _, p := range dp {
		for _, g := range dg {
			pgfn2 += (p * g) * (p * g)
		}
	}
	egfn2 := 0.0
	for _, e := range de {
		for _, g := range dg {
			egfn2 += (e * g) * (e * g)
		}
	}

	multicell.ProjectSVD("<dp_dG>", dirE, dirG, deltaP, deltaG, false, false)
	fmt.Printf("<dp><dG>_FN2,Tr\t\t%e\t%e\n", pgfn2, 0)
	multicell.ProjectSVD("<Ddp_DdG>", dirE, dirG, deltaP, deltaG, true, true)
	multicell.ProjectSVD("  <Dp1_DG1>", dirE, dirG, p1, G1, true, true)
	multicell.ProjectSVD("  <Dp0_DG0>", dirE, dirG, p0, G0, true, true)
	multicell.ProjectSVD("  <Dp0_DG1>", dirE, dirG, p0, G1, true, true)
	multicell.ProjectSVD("  <Dp1_DG0>", dirE, dirG, p1, G0, true, true)

	multicell.ProjectSVD("<de_dG>", dirE, dirG, deltaE, deltaG, false, false)
	fmt.Printf("<de><dG>_FN2,Tr\t\t%e\t%e\n", egfn2, 0)
	multicell.ProjectSVD("<Dde_DdG>", dirE, dirG, deltaE, deltaG, true, true)
	multicell.ProjectSVD("  <De1_DG1>", dirE, dirG, e1, G1, true, true)
	multicell.ProjectSVD("  <De0_DG0>", dirE, dirG, e0, G0, true, true)
	multicell.ProjectSVD("  <De0_DG1>", dirE, dirG, e0, G1, true, true)
	multicell.ProjectSVD("  <De1_DG0>", dirE, dirG, e1, G0, true, true)

}
