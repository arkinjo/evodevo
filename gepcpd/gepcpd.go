package main

// Don't use this: it is not very useful.
// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/arkinjo/evodevo/multicell"
	//	"gonum.org/v1/gonum/mat"
)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	jsonP := flag.String("jsonin", "", "json file of population")
	maxiterP := flag.Int("maxiter", 100, "Maximum number of CPD iterations")
	rankP := flag.Int("rank", 41, "Maximum number of CPD iterations")
	flag.Parse()

	maxiter := *maxiterP
	rank := *rankP

	pop := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsonP != "" {
		pop.ImportPopGz(*jsonP)
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

	//	mg, me, mp, cov := multicell.GetCrossCov3(genome, dele, delp, true, true, true)
	mg, me, mp, cov := MakeFakeCov3(genome, dele, delp, rank)

	dirg := multicell.CopyVec(mg)
	dire := multicell.CopyVec(me)
	dirp := multicell.CopyVec(mp)
	multicell.ScaleVec(dirg, 1.0/multicell.Norm2(dirg), dirg)
	multicell.ScaleVec(dire, 1.0/multicell.Norm2(dire), dire)
	multicell.ScaleVec(dirp, 1.0/multicell.Norm2(dirp), dirp)

	cpd := multicell.GetCPDO(cov, denv, denv, maxiter, rank)

	dev, rat := multicell.CheckCPD(cov, cpd)
	fmt.Printf("CPD_devrat\t%d\t%e\t%e\n", len(cpd), dev, rat)

	for r, p := range cpd {
		fmt.Printf("CPD_vals\t%d\t%e\n", r, p.SVal)
	}

	for a := range cpd {
		dotg := multicell.DotVecs(dirg, cpd[a].Axes[0])
		dote := multicell.DotVecs(dire, cpd[a].Axes[1])
		dotp := multicell.DotVecs(dirp, cpd[a].Axes[2])
		fmt.Printf("CPD_ali\t%d\t%e\t%e\t%e\n", a, dotg, dote, dotp)
	}

	ncol := rank
	if ncol > 4 {
		ncol = 4
	}
	for k, g := range genome {
		fmt.Printf("CPD_prj\t%d", k)
		e := dele[k]
		p := delp[k]
		multicell.DiffVecs(g, g, mg)
		multicell.DiffVecs(e, e, me)
		multicell.DiffVecs(p, p, mp)

		for a := 0; a < ncol; a++ {
			dg := multicell.DotVecs(g, cpd[a].Axes[0])
			de := multicell.DotVecs(e, cpd[a].Axes[1])
			dp := multicell.DotVecs(p, cpd[a].Axes[2])
			fmt.Printf("\t%e\t%e\t%e", dg, de, dp)

		}
		fmt.Println("")
	}
}

func MakeFakeCov3(genome, dele, delp multicell.Dmat, rank int) ([]float64, []float64, []float64, multicell.Tensor3) {
	mg, me, cov0 := multicell.GetCrossCov(genome, dele, true, true)
	u0, s0, u1 := multicell.GetSVD(cov0)

	_, mp, cov1 := multicell.GetCrossCov(dele, delp, true, true)
	_, s1, u2 := multicell.GetSVD(cov1)

	len0, len1, len2 := len(mg), len(me), len(mp)
	ten := multicell.NewTensor3(len0, len1, len2)

	for i := 0; i < len0; i++ {
		for j := 0; j < len1; j++ {
			for k := 0; k < len2; k++ {
				for r := 0; r < rank; r++ {
					ten[i][j][k] += math.Sqrt(s0[r]*s1[r]) * u0.At(i, r) * u1.At(j, r) * u2.At(k, r)
				}
			}
		}
	}
	return mg, me, mp, ten
}
