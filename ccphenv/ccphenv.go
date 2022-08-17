package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"

	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")

	jsongzinP := flag.String("jsonin", "", ".json.gz file of input population")

	flag.Parse()

	pop := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsongzinP != "" {
		pop.ImportPopGz(*jsongzinP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	dim := len(env0)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)

	e0 := pop.GetFlatStateVec("E", 0)
	e1 := pop.GetFlatStateVec("E", 1)
	me0 := multicell.GetMeanVec(e0)
	me1 := multicell.GetMeanVec(e1)
	de := multicell.NewVec(dim)
	multicell.DiffVecs(de, me1, me0)

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)
	mp0 := multicell.GetMeanVec(p0)
	mp1 := multicell.GetMeanVec(p1)
	dp := multicell.NewVec(dim)
	multicell.DiffVecs(dp, mp1, mp0)

	deltaP := make([]multicell.Vec, 0)
	deltaE := make([]multicell.Vec, 0)
	for k, p := range p1 {
		d := multicell.NewVec(dim)
		multicell.DiffVecs(d, p, p0[k])
		deltaP = append(deltaP, d)

		e := multicell.NewVec(dim)
		multicell.DiffVecs(e, e1[k], e0[k])
		deltaE = append(deltaE, e)
	}

	dotPE := multicell.DotVecs(dp, de)
	dotPP := multicell.DotVecs(dp, dp)
	dotEE := multicell.DotVecs(de, de)

	dirE := mat.NewVecDense(dim, multicell.CopyVec(denv))
	dirE.ScaleVec(1.0/dirE.Norm(2), dirE)
	dirP := mat.NewVecDense(dim, multicell.CopyVec(dp))
	dirP.ScaleVec(1.0/dirP.Norm(2), dirP)

	Project("<dp_dp>", dirE, dirP, deltaP, deltaP, false, false)
	fmt.Printf("||<dp>||^2\t%e\n", dotPP)
	Project("<Ddp_Ddp>", dirE, dirP, deltaP, deltaP, true, true)
	Project(" <Dp1_Dp1>", dirE, dirP, p1, p1, true, true)
	Project(" <Dp0_Dp0>", dirE, dirP, p0, p0, true, true)
	Project(" <Dp1_Dp0>", dirE, dirP, p1, p0, true, true)
	Project(" <Dp0_Dp1>", dirE, dirP, p0, p1, true, true)

	Project("<dp_de>", dirE, dirP, deltaP, deltaE, false, false)
	fmt.Printf("<dp><de>FN2,Tr\t\t%e\t%e\n", dotPP*dotEE, dotPE)
	Project("<Ddp_Dde>", dirE, dirP, deltaP, deltaE, true, true)
	Project(" <Dp1_De1>", dirE, dirP, p1, e1, true, true)
	Project(" <Dp0_De0>", dirE, dirP, p0, e0, true, true)
	Project(" <Dp1_De0>", dirE, dirP, p1, e0, true, true)
	Project(" <Dp0_De1>", dirE, dirP, p0, e1, true, true)

	Project("<de_de>", dirE, dirE, deltaE, deltaE, false, false)
	fmt.Printf("||<de>||^2\t%e\n", dotEE)
	Project("<Dde_Dde>", dirE, dirP, deltaE, deltaE, true, true)
	Project(" <De1_De1>", dirE, dirP, e1, e1, true, true)
	Project(" <De0_De0>", dirE, dirP, e0, e0, true, true)
	Project(" <De1_De0>", dirE, dirP, e1, e0, true, true)
	Project(" <De0_De1>", dirE, dirP, e0, e1, true, true)

	mixp := make([]multicell.Vec, 0)
	mixe := make([]multicell.Vec, 0)
	for k := range p0 {
		mixp = append(mixp, p0[k])
		mixe = append(mixe, e0[k])
	}
	for k := range p1 {
		mixp = append(mixp, p1[k])
		mixe = append(mixe, e1[k])
	}
	Project("mixPP", dirE, dirP, mixp, mixp, true, true)
	Project("mixPE", dirE, dirP, mixp, mixe, true, true)

	LinearResponse(&pop)
}

func Project(label string, dir0, dir1 *mat.VecDense, data0, data1 [][]float64, sub0, sub1 bool) {
	mean0, mean1, ccmat := multicell.GetCrossCov(data0, data1, sub0, sub1)
	U, vals, V := multicell.GetSVD(ccmat)
	multicell.ProjectSVD(label, dir0, dir1, data0, data1, mean0, mean1, ccmat, vals, U, V)
}

func LinearResponse(pop *multicell.Population) {
	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	dim := len(env0)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)

	fe0 := pop.GetFlatStateVec("E", 0)
	fp0 := pop.GetFlatStateVec("P", 0)
	fp1 := pop.GetFlatStateVec("P", 1)
	p0_ave, _, pe_cov := multicell.GetCrossCov(fp0, fe0, true, true)
	p1_ave := multicell.GetMeanVec(fp1)
	dp := multicell.NewVec(dim)
	multicell.DiffVecs(dp, p1_ave, p0_ave)

	dpp := multicell.NewVec(dim)
	sigma := pop.Params.SDNoise
	sigma2 := (sigma * sigma)
	for i, ci := range pe_cov {
		for j, v := range ci {
			dpp[i] += v * denv[j]
		}
	}

	for i, de := range denv {
		fmt.Printf("LRT\t%2d\t%e\t%e\t%e\n", i, de, dp[i], dpp[i]/sigma2)
	}
}
