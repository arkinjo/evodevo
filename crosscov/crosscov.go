package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/arkinjo/evodevo/multicell"
)

var sqrt3 float64 = math.Sqrt(3.0)

func cca(vecs0, vecs1 [][]float64, axis0 []float64) {
	len0 := len(vecs0[0])
	len1 := len(vecs1[0])

	av0, _, cov := multicell.GetCrossCov(vecs0, vecs1, true, true)
	Us, svals, Vs := multicell.GetSVD(cov)

	rank := len(svals)
	U := multicell.NewDmat(rank, len0)
	V := multicell.NewDmat(rank, len1)

	for a := range svals {
		for i := range U[a] {
			U[a][i] = Us.At(i, a)
		}
		for i := range V[a] {
			V[a][i] = Vs.At(i, a)
		}

		neg := multicell.DotVecs(U[a], axis0)
		if neg < 0.0 {
			multicell.ScaleVec(U[a], -1, U[a])
			multicell.ScaleVec(V[a], -1, V[a])
		}
	}

	totS := 0.0
	for _, v := range svals {
		totS += v * v
	}
	fmt.Printf("FN2\t%e\n", totS)
	cum := 0.0
	for a, v := range svals {
		v2 := v * v
		cum += v * v / totS
		fmt.Printf("SVals\t%d\t%e\t%e\t%e\n", a, v2, v2/totS, cum) //ath singular value, prop, cum prop
		ali := multicell.DotVecs(U[a], axis0)
		fmt.Printf("Ali\t%d\t%e\n", a, ali)
	}
	for i, v := range U[0] {
		fmt.Printf("LSVec1\t%d\t%e\n", i, v)
	}

	// Projecting individuals' phenotypes onto U[0].
	dp := multicell.NewVec(len0)
	for i, p := range vecs0 {
		multicell.DiffVecs(dp, p, av0)
		v := multicell.DotVecs(dp, U[0])
		fmt.Printf("LProj1\t%d\t%e\n", i, v)
	}
}

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	jsonP := flag.String("jsonin", "", "json file of population")
	modeP := flag.Int("mode", 2, "0: AncEnv; 1: NovEnv; 2: NovEnv - AncEnv")
	rflagP := flag.String("right", "e", "Right component (e or G)")
	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	modeFlag := *modeP
	rightFlag := *rflagP

	switch modeFlag {
	case 0:
		fmt.Println("# Ancestral Environment")
	case 1:
		fmt.Println("# Novel Environment")
	case 2:
		fmt.Println("# difference between Novel and Ancestral Environments")
	default:
		log.Fatal("-mode must be 0, 1, or 2")
	}

	pop := multicell.NewPopulation(settings)
	if *jsonP != "" {
		pop.ImportPopGz(*jsonP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	Nenv := multicell.GetNenv()
	Nsel := multicell.GetNsel()

	fenv0 := multicell.FlattenEnvs(pop.AncEnvs)
	//	fenv1 := multicell.FlattenEnvs(pop.NovEnvs)
	lenE := len(fenv0)

	env0 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop.AncEnvs))
	env1 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop.NovEnvs))

	lenP := len(env0)
	denv := multicell.NewVec(lenP)
	multicell.DiffVecs(denv, env1, env0)
	// for PC projection
	paxis := multicell.CopyVec(denv) //paxis is defined as change in environment
	multicell.NormalizeVec(paxis)

	genome0 := pop.GetFlatGenome(multicell.IAncEnv)
	genome1 := pop.GetFlatGenome(multicell.INovEnv)
	lenG := len(genome0[0])
	delg := make([][]float64, 0)
	for k, g := range genome0 {
		switch modeFlag {
		case 0:
			delg = append(delg, g)
		case 1:
			delg = append(delg, genome1[k])
		case 2:
			d := multicell.NewVec(lenG)
			multicell.DiffVecs(d, genome1[k], g)
			delg = append(delg, d)
		}

	}

	e0 := pop.GetFlatStateVec("E", 0, 0, Nenv)
	e1 := pop.GetFlatStateVec("E", 1, 0, Nenv)
	dele := make([][]float64, 0)
	for k, e := range e0 {
		switch modeFlag {
		case 0:
			dele = append(dele, e)
		case 1:
			dele = append(dele, e1[k])
		case 2:
			d := multicell.NewVec(lenE)
			multicell.DiffVecs(d, e1[k], e)
			dele = append(dele, d)
		}

	}

	p0 := pop.GetFlatStateVec("P", 0, 0, Nsel) //ancestral phenotype
	p1 := pop.GetFlatStateVec("P", 1, 0, Nsel) //novel phenotype
	delp := make([][]float64, 0)
	for k, p := range p0 {
		switch modeFlag {
		case 0:
			delp = append(delp, p)
		case 1:
			delp = append(delp, p1[k])
		case 2:
			d := multicell.NewVec(lenP)
			multicell.DiffVecs(d, p1[k], p)
			delp = append(delp, d)
		}

	}

	if rightFlag == "e" {
		cca(delp, dele, paxis)
	} else if rightFlag == "G" {
		cca(delp, delg, paxis)
	} else {
		log.Fatal("-right must be either e or G.")
	}
}
