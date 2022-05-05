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
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")

	jsonP := flag.String("jsonin", "", "json file of population")
	jsonrefP := flag.String("jsonref", "", "json file of reference population")

	flag.Parse()

	pop := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsonP != "" {
		pop.FromJSON(*jsonP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsonrefP != "" {
		pop1.FromJSON(*jsonrefP)
		multicell.SetParams(pop1.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the reference JSON file with -jsonref=filename.")
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	lenE := len(env0)
	denv := multicell.NewVec(lenE)
	multicell.DiffVecs(denv, env1, env0)
	multicell.ScaleVec(denv, 1.0/multicell.Norm2(denv), denv)

	genome := pop.GetFlatGenome()
	genome1 := pop1.GetFlatGenome()
	gref := multicell.GetMeanVec(genome1)
	lenG := len(genome[0])
	delg := make([][]float64, 0)
	for _, g := range genome {
		d := multicell.NewVec(lenG)
		multicell.DiffVecs(d, g, gref)
		delg = append(delg, d)
	}
	e0 := pop.GetFlatStateVec("E", 0)
	e1 := pop.GetFlatStateVec("E", 1)
	dele := make([][]float64, 0)
	for k, e := range e0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, e1[k], e)
		dele = append(dele, d)
	}

	deleg := make([][]float64, 0)
	for k, e := range dele {
		d := append(e, delg[k]...)
		deleg = append(deleg, d)
	}

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)
	delp := make([][]float64, 0)
	for i, p := range p0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, p1[i], p)
		delp = append(delp, d)
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
	Up := multicell.NewDmat(rank, lenE)
	Veg := multicell.NewDmat(rank, lenE+lenG)
	Ve := multicell.NewDmat(rank, lenE)
	Vg := multicell.NewDmat(rank, lenG)
	weightE := multicell.NewVec(rank)

	for a := range svals {
		k := 0
		mage := 0.0
		magg := 0.0
		for i := 0; i < lenE; i++ {
			Up[a][i] = U.At(i, a)
			v := V.At(i, a)
			Ve[a][i] = v
			Veg[a][i] = v
			mage += v * v
			k++
		}
		for i := 0; i < lenG; i++ {
			v := V.At(i, a)
			Vg[a][i] = v
			Veg[a][lenE+i] = v
			magg += v * v
			k++
		}
		weightE[a] = mage
		multicell.NormalizeVec(Ve[a])
		multicell.NormalizeVec(Vg[a])

		neg := multicell.DotVecs(Up[a], denv)
		if neg < 0.0 {
			multicell.ScaleVec(Up[a], -1, Up[a])
			multicell.ScaleVec(Ve[a], -1, Ve[a])
			multicell.ScaleVec(Vg[a], -1, Vg[a])
			multicell.ScaleVec(Veg[a], -1, Veg[a])
		}
	}

	fne := 0.0
	fng := 0.0
	for _, p := range mp {
		for _, e := range me {
			fne += (p * e) * (p * e)
		}
		for _, g := range mg {
			fng += (p * g) * (p * g)
		}
	}
	fmt.Printf("<dp><de>_FN2\t%e\n", fne)
	fmt.Printf("<dp><dG>_FN2\t%e\n", fng)
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
		dotg := multicell.DotVecs(mg, Vg[a])
		dote := multicell.DotVecs(me, Ve[a])
		dotp := multicell.DotVecs(mp, Up[a])

		fmt.Printf("\t%e\t%e\t%e\t%e", doteg, dotg, dote, dotp)
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

	for i := 0; i < lenE; i++ {
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
