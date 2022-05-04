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

	mg, me, mp, cov := multicell.GetCrossCov3(genome, dele, delp, true, true, true)
	dirg := multicell.CopyVec(mg)
	dire := multicell.CopyVec(me)
	dirp := multicell.CopyVec(mp)
	multicell.ScaleVec(dirg, 1.0/multicell.Norm2(dirg), dirg)
	multicell.ScaleVec(dire, 1.0/multicell.Norm2(dire), dire)
	multicell.ScaleVec(dirp, 1.0/multicell.Norm2(dirp), dirp)

	cpd := multicell.GetCPDO(cov, 100)

	for i, p := range cpd {
		fmt.Printf("CPD_vals %d %d %e\n", i, p.I, p.SVal)
	}

	for a := range cpd {
		dotg := multicell.DotVecs(dirg, cpd[a].Axes[0])
		dote := multicell.DotVecs(dire, cpd[a].Axes[1])
		dotp := multicell.DotVecs(dirp, cpd[a].Axes[2])
		fmt.Printf("CPD_ali\t%d\t%e\t%e\t%e\n", a, dotg, dote, dotp)
	}

	for k, g := range genome {
		fmt.Printf("CPD_prj\t%d", k)
		e := dele[k]
		p := delp[k]
		multicell.DiffVecs(g, g, mg)
		multicell.DiffVecs(e, e, me)
		multicell.DiffVecs(p, p, mp)
		for a := 0; a < 3; a++ {
			dg := multicell.DotVecs(g, cpd[a].Axes[0])
			de := multicell.DotVecs(e, cpd[a].Axes[1])
			dp := multicell.DotVecs(p, cpd[a].Axes[2])
			fmt.Printf("\t%e\t%e\t%e", dg, de, dp)
		}
		fmt.Println("")
	}
}
