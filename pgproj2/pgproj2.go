package main

import (
	"flag"
	"fmt"
	"github.com/arkinjo/evodevo/multicell"
	"log"
	"os"
	"time"
)

var PG_Filename string //Dump for phenotypes and genotypes
var json_in string     //JSON encoding of initial population; default to empty string

func main() {
	log.Println("Starting...")
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	ref1Ptr := flag.String("ref1", "", "reference JSON file 1")
	ref2Ptr := flag.String("ref2", "", "reference JSON file 2")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "basename of JSON files") //default to empty string

	flag.Parse()

	epochlength := *genPtr
	refgen1 := *ref1Ptr
	refgen2 := *ref2Ptr
	PG_Filename = *pgfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	if refgen1 == "" {
		log.Fatal("Must specify JSON reference file 1.")
	}

	log.Println("epochlength", epochlength)

	tdump := time.Now()

	log.Println("Reading Pop1")
	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP)
	fmt.Println("Reference population :", refgen1)
	pop1.FromJSON(refgen1)
	multicell.SetParams(pop1.Params)
	// Reference direction
	env0 := multicell.FlattenEnvs(pop1.AncEnvs)
	env1 := multicell.FlattenEnvs(pop1.NovEnvs)
	lenE := len(env0)
	denv := multicell.NewVec(lenE)
	multicell.DiffVecs(denv, env1, env0)
	multicell.NormalizeVec(denv)

	g1 := pop1.GetFlatGenome()
	e10 := pop1.GetFlatStateVec("E", 0)
	e11 := pop1.GetFlatStateVec("E", 1)
	for k, g := range g1 {
		e10[k] = append(e10[k], g...)
		e11[k] = append(e11[k], g...)
	}
	p10 := pop1.GetFlatStateVec("P", 0)
	p11 := pop1.GetFlatStateVec("P", 1)

	gmix := make([][]float64, 0)
	pmix := make([][]float64, 0)

	gmix = append(gmix, e10...)
	pmix = append(pmix, p10...)
	gmix = append(gmix, e11...)
	pmix = append(pmix, p11...)
	log.Println("lenG,lenP=", len(gmix[0]), len(pmix[0]))
	if refgen2 != "" {
		log.Println("Reading Pop2")
		pop2 := multicell.NewPopulation(*ncellsP, *maxpopP)
		fmt.Println("Reference population 2:", refgen2)
		pop2.FromJSON(refgen2)

		g2 := pop2.GetFlatGenome()
		e20 := pop2.GetFlatStateVec("E", 0)
		e21 := pop2.GetFlatStateVec("E", 1)
		for k, g := range g2 {
			e20[k] = append(e20[k], g...)
			e21[k] = append(e21[k], g...)
		}
		p20 := pop2.GetFlatStateVec("P", 0)
		p21 := pop2.GetFlatStateVec("P", 1)
		gmix = append(gmix, e20...)
		pmix = append(pmix, p20...)
		gmix = append(gmix, e21...)
		pmix = append(pmix, p21...)
	}

	log.Println("Finding Principal Axes")
	mp, mg, cov := multicell.GetCrossCov(pmix, gmix, true, true)
	U, sval, V := multicell.GetSVD(cov)
	log.Println("Svals:", sval)
	paxis := multicell.NewDmat(2, len(mp))
	gaxis := multicell.NewDmat(2, len(mg))
	for a, pa := range paxis {
		for i := range pa {
			paxis[a][i] = U.At(i, a)
		}
	}
	for a, ga := range gaxis {
		for i := range ga {
			gaxis[a][i] = V.At(i, a)
		}
	}

	for a, pa := range paxis {
		neg := multicell.DotVecs(pa, denv)
		if neg < 0.0 {
			multicell.ScaleVec(paxis[a], -1.0, pa)
			multicell.ScaleVec(gaxis[a], -1.0, gaxis[a])
		}
	}
	log.Printf("Dumping start")
	for gen := 1; gen <= epochlength; gen++ {
		ofilename := fmt.Sprintf("%s_%3.3d.dat", PG_Filename, gen)
		fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintf(fout, "#Geno+e0(0)  \tPheno0(0)  \tGeno+e1(0)  \tPheno1(0)")
		fmt.Fprintf(fout, "\tGeno+e0(1)  \tPheno0(1)  \tGeno+e1(1)  \tPheno1(1)")
		fmt.Fprintf(fout, "\t||p0-e0||  \t||p1-e1||  \tFit     \tWagFit\n")

		jfilename := fmt.Sprintf("%s_%3.3d.json", json_in, gen)
		pop := multicell.NewPopulation(*ncellsP, *maxpopP)
		pop.FromJSON(jfilename)
		gt := pop.GetFlatGenome()
		et0 := pop.GetFlatStateVec("E", 0)
		et1 := pop.GetFlatStateVec("E", 1)
		for k, g := range gt {
			et0[k] = append(et0[k], g...)
			et1[k] = append(et1[k], g...)
		}

		pt0 := pop.GetFlatStateVec("P", 0)
		pt1 := pop.GetFlatStateVec("P", 1)

		for k := range et0 {
			multicell.DiffVecs(et0[k], et0[k], mg)
			multicell.DiffVecs(et1[k], et1[k], mg)
			multicell.DiffVecs(pt0[k], pt0[k], mp)
			multicell.DiffVecs(pt1[k], pt1[k], mp)

			for a, pa := range paxis {
				x0 := multicell.DotVecs(et0[k], gaxis[a])
				y0 := multicell.DotVecs(pt0[k], pa)
				x1 := multicell.DotVecs(et1[k], gaxis[a])
				y1 := multicell.DotVecs(pt1[k], pa)
				fmt.Fprintf(fout, "\t%e\t%e\t%e\t%e", x0, y0, x1, y1)
			}
			dp1e1 := pop.Indivs[k].Dp1e1
			dp0e0 := pop.Indivs[k].Dp0e0
			fit := pop.Indivs[k].Fit
			wf := pop.Indivs[k].WagFit

			fmt.Fprintf(fout, "\t%e\t%e\t%e\t%e\n", dp0e0, dp1e1, fit, wf)

		}
		err = fout.Close()
		if err != nil {
			log.Fatal(err)
		}

	}

	dtdump := time.Since(tdump)
	fmt.Println("Time taken to dump projections :", dtdump)
	fmt.Printf("Projections written to %s_*.dat \n", PG_Filename)
	dt := time.Since(t0)

	fmt.Println("Total time taken : ", dt)
}
