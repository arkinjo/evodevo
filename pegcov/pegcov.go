package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/arkinjo/evodevo/multicell"
	//	"gonum.org/v1/gonum/mat"
)

var sqrt3 float64 = math.Sqrt(3.0)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	jsonP := flag.String("jsonin", "", "json file of population")
	pflagP := flag.Int("p", 2, "0: AncEnv; 1: NovEnv; 2: NovEnv - AncEnv")
	egflagP := flag.Int("eg", 2, "0: AncEnv; 1: NovEnv; 2: NovEnv - AncEnv")
	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	pFlag := *pflagP
	egFlag := *egflagP

	switch pFlag {
	case 0:
		fmt.Println("# phenotype under Ancestral Environment")
	case 1:
		fmt.Println("# phenotype under Novel Environment")
	case 2:
		fmt.Println("# phenotype difference between Novel and Ancestral Environments")
	default:
		log.Fatal("-p must be 0, 1, or 2")
	}
	switch egFlag {
	case 0:
		fmt.Println("# genotype (+env) under Ancestral Environment")
	case 1:
		fmt.Println("# genotype (+env) under Novel Environment")
	case 2:
		fmt.Println("# genotype (+env) difference between Novel and Ancestral Environments")
	default:
		log.Fatal("-p must be 0, 1, or 2")
	}

	pop := multicell.NewPopulation(settings)
	if *jsonP != "" {
		pop.FromJSON(*jsonP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	Nenv := multicell.GetNenv()
	Nsel := multicell.GetNsel()

	env0 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop.AncEnvs))
	env1 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop.NovEnvs))
	lenE := len(env0)
	lenP := multicell.GetNsel()
	denv := multicell.NewVec(lenP)
	multicell.DiffVecs(denv, env1[0:lenP], env0[0:lenP])
	multicell.NormalizeVec(denv)

	genome0 := pop.GetFlatGenome(multicell.IAncEnv)
	genome1 := pop.GetFlatGenome(multicell.INovEnv)
	lenG := len(genome0[0])
	delg := make([][]float64, 0)
	for k, g := range genome0 {
		switch egFlag {
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
		switch egFlag {
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

	deleg := make([][]float64, 0)
	for k, e := range dele {
		d := append(e, delg[k]...)
		deleg = append(deleg, d)
	}

	p0 := pop.GetFlatStateVec("P", 0, 0, Nsel)
	p1 := pop.GetFlatStateVec("P", 1, 0, Nsel)
	delp := make([][]float64, 0)
	for k, p := range p0 {
		switch pFlag {
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

	mp, meg, cov := multicell.GetCrossCov(delp, deleg, true, true)

	me := meg[0:lenE]
	mg := meg[lenE:]
	for k, p := range delp {
		multicell.DiffVecs(delp[k], p, mp)
		multicell.DiffVecs(dele[k], dele[k], me)
		multicell.DiffVecs(delg[k], delg[k], mg)
		multicell.DiffVecs(deleg[k], deleg[k], meg)
	}
	U, svals, V := multicell.GetSVD(cov)
	rank := len(svals)
	Up := multicell.NewDmat(rank, lenP)
	Veg := multicell.NewDmat(rank, lenE+lenG)
	Ve := multicell.NewDmat(rank, 0)
	Vg := multicell.NewDmat(rank, 0)
	weightE := multicell.NewVec(rank)

	for a := range svals {
		for i := range Up[a] {
			Up[a][i] = U.At(i, a)
		}
		for i := range Veg[a] {
			Veg[a][i] = V.At(i, a)
		}

		neg := multicell.DotVecs(Up[a], denv)
		if neg < 0.0 {
			multicell.ScaleVec(Up[a], -1, Up[a])
			multicell.ScaleVec(Veg[a], -1, Veg[a])
		}
		Ve[a] = multicell.CopyVec(Veg[a][0:lenE])
		Vg[a] = multicell.CopyVec(Veg[a][lenE:])
		weightE[a] = multicell.Norm2Sq(Ve[a])
		multicell.NormalizeVec(Ve[a])
		multicell.NormalizeVec(Vg[a])
	}

	fne := 0.0
	fng := 0.0
	fnp := 0.0
	for _, p := range mp {
		fnp += p * p
	}
	for _, e := range me {
		fne += e * e
	}
	for _, g := range mg {
		fng += g * g
	}

	fmt.Printf("<dp><de><dG>_FN2\t%e\t%e\t%e\n", fnp, fne, fng)

	totS := 0.0
	totSe := 0.0
	for a, v := range svals {
		totS += v * v
		totSe += v * v * weightE[a]
	}
	fmt.Printf("<DdpDde>_FN2\t%e\n", totSe)
	fmt.Printf("<DdpDdG>_FN2\t%e\n", totS-totSe)
	cum := 0.0
	for a, v := range svals {
		v2 := v * v
		cum += v * v / totS
		fmt.Printf("SVals\t%d\t%e\t%e\t%e\t%e\n", a, v2, v2/totS, cum, weightE[a])
		ali := multicell.DotVecs(Up[a], denv)
		fmt.Printf("Ali\t%d\t%e\n", a, ali)
	}
	fmt.Printf("Center\t0")
	for a := 0; a < 3; a++ {
		doteg := multicell.DotVecs(meg, Veg[a])
		dote := multicell.DotVecs(me, Ve[a])
		dotg := multicell.DotVecs(mg, Vg[a])
		dotp := multicell.DotVecs(mp, Up[a])

		fmt.Printf("\t%e\t%e\t%e\t%e", doteg, dote, dotg, dotp)
	}
	fmt.Println()

	for k, dp0 := range delp {
		fmt.Printf("Prj\t%d", k)
		dg0 := delg[k]
		de0 := dele[k]
		deg0 := deleg[k]
		for i := 0; i < 3; i++ {
			doteg := multicell.DotVecs(deg0, Veg[i])
			dotg := multicell.DotVecs(dg0, Vg[i])
			dote := multicell.DotVecs(de0, Ve[i])
			dotp := multicell.DotVecs(dp0, Up[i])

			fmt.Printf("\t%e\t%e\t%e\t%e", doteg, dote, dotg, dotp)
		}
		fmt.Println()
	}

	for i := 0; i < lenP; i++ {
		fmt.Printf("<dp>,Uvec\t%d\t%e", i, mp[i])
		for a := 0; a < 5; a++ {
			fmt.Printf("\t%e", Up[a][i])
		}
		fmt.Println()
	}
	for i := 0; i < lenE+lenG; i++ {
		fmt.Printf("<de+dG>,Vvec\t%d\t%e", i, meg[i])
		for a := 0; a < 5; a++ {
			fmt.Printf("\t%e", Veg[a][i])
		}
		fmt.Println()
	}
}
