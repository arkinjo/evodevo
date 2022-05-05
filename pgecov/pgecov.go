package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
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
	multicell.ScaleVec(denv, 1.0/multicell.Norm2(denv), denv)
	dirE := mat.NewVecDense(lenE, denv)

	genome := pop.GetFlatGenome()
	mg := multicell.GetMeanVec(genome)
	lenG := len(genome[0])
	delg := make([][]float64, 0)
	for _, g := range genome {
		d := multicell.NewVec(lenG)
		multicell.DiffVecs(d, g, mg)
		delg = append(delg, d)
	}

	e0 := pop.GetFlatStateVec("E", 0)
	//	me0 := multicell.GetMeanVec(e0)
	e1 := pop.GetFlatStateVec("E", 1)
	dele := make([][]float64, 0)
	for k, e := range e0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, e1[k], e)
		//	multicell.DiffVecs(d, d, me0)
		dele = append(dele, d)
	}
	//	me := multicell.GetMeanVec(dele)
	deleg := make([][]float64, 0)
	for k, e := range dele {
		d := append(e, delg[k]...)
		deleg = append(deleg, d)
	}
	//	meg := multicell.GetMeanVec(deleg)

	p0 := pop.GetFlatStateVec("P", 0)
	//	mp0 := multicell.GetMeanVec(p0)
	p1 := pop.GetFlatStateVec("P", 1)
	delp := make([][]float64, 0)
	for i, p := range p0 {
		d := multicell.NewVec(lenE)
		multicell.DiffVecs(d, p1[i], p)
		//	multicell.DiffVecs(d, d, mp0)
		delp = append(delp, d)
	}
	//	mp := multicell.GetMeanVec(delp)
	_, _, cov := multicell.GetCrossCov(delp, deleg, false, false)
	U, svals, V := multicell.GetSVD(cov)
	rank := len(svals)
	Up := multicell.NewDmat(rank, lenE)
	Ve := multicell.NewDmat(rank, lenE)
	Vg := multicell.NewDmat(rank, lenG)
	weightE := multicell.NewVec(rank)

	for a := range svals {
		k := 0
		mage := 0.0
		magg := 0.0
		for i := 0; i < lenE; i++ {
			v := V.At(i, a)
			Up[a][i] = U.At(i, a)
			Ve[a][i] = v
			mage += v * v
			k++
		}
		for i := 0; i < lenG; i++ {
			v := V.At(i, a)
			Vg[a][i] = v
			magg += v * v
			k++
		}
		weightE[a] = mage
		scae := 1.0 / math.Sqrt(mage)
		scag := 1.0 / math.Sqrt(magg)
		multicell.ScaleVec(Ve[a], scae, Ve[a])
		multicell.ScaleVec(Vg[a], scag, Vg[a])

		neg := multicell.DotVecs(Up[a], denv)
		if neg < 0.0 {
			for i := 0; i < lenE; i++ {
				Up[a][i] *= -1.0
				Ve[a][i] *= -1.0
			}
			for i := 0; i < lenG; i++ {
				Vg[a][i] *= -1.0
			}
		}
	}

	totS := 0.0
	for _, v := range svals {
		totS += v
	}
	cum := 0.0
	for i, v := range svals {
		cum += v / totS
		fmt.Printf("SVals\t%d\t%e\t%e\t%e\t%e\n", i, v, v/totS, cum, weightE[i])
		ui := U.ColView(i)
		ali := mat.Dot(ui, dirE)
		fmt.Printf("Ali\t%d\t%e\n", i, ali)
	}
	for k, dp0 := range delp {
		fmt.Printf("Prj\t%d", k)
		dg0 := delg[k]
		de0 := dele[k]

		for i := 0; i < 3; i++ {
			dotg := multicell.DotVecs(dg0, Vg[i])
			dote := multicell.DotVecs(de0, Ve[i])
			dotp := multicell.DotVecs(dp0, Up[i])

			fmt.Printf("\t%e\t%e\t%e", dotg, dote, dotp)
		}
		fmt.Println()
	}
}
