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
	ref1Ptr := flag.String("ref1", "", "reference JSON file 1 (requried)")
	ref2Ptr := flag.String("ref2", "", "reference JSON file 2 (optional)")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "basename of JSON files")

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
	if refgen2 == "" {
		log.Fatal("Must specify JSON reference file 2.")
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

	g10 := pop1.GetFlatGenome(multicell.IAncEnv)
	e10 := pop1.GetFlatStateVec("E", 0)
	for k, g := range g10 {
		e10[k] = append(e10[k], g...)
	}

	log.Println("Reading Pop2")
	pop2 := multicell.NewPopulation(*ncellsP, *maxpopP)
	fmt.Println("Reference population 2:", refgen2)
	pop2.FromJSON(refgen2)

	g21 := pop2.GetFlatGenome(multicell.INovEnv)
	e21 := pop2.GetFlatStateVec("E", 1)
	for k, g := range g21 {
		e21[k] = append(e21[k], g...)
	}

	log.Println("Finding Principal Axes")
	mp := multicell.NewVec(lenE)
	multicell.AddVecs(mp, env0, env1)
	multicell.ScaleVec(mp, 0.5, mp)
	paxis := multicell.CopyVec(denv)

	mg0 := multicell.GetMeanVec(e10)
	mg1 := multicell.GetMeanVec(e21)
	lenG := len(mg0)
	mg := multicell.NewVec(lenG)
	multicell.AddVecs(mg, mg0, mg1)
	multicell.ScaleVec(mg, 0.5, mg)

	gaxis := multicell.NewVec(lenG)
	multicell.DiffVecs(gaxis, mg1, mg0)
	multicell.NormalizeVec(gaxis)

	{
		et0 := multicell.NewVec(lenG)
		et1 := multicell.NewVec(lenG)
		multicell.DiffVecs(et0, mg0, mg)
		multicell.DiffVecs(et1, mg1, mg)
		pt0 := multicell.NewVec(lenE)
		pt1 := multicell.NewVec(lenE)
		multicell.DiffVecs(pt0, env0, mp)
		multicell.DiffVecs(pt1, env1, mp)
		x0 := multicell.DotVecs(et0, gaxis) / multicell.Norm2(et0)
		y0 := multicell.DotVecs(pt0, paxis) / multicell.Norm2(pt0)
		x1 := multicell.DotVecs(et1, gaxis) / multicell.Norm2(et1)
		y1 := multicell.DotVecs(pt1, paxis) / multicell.Norm2(pt1)
		log.Println("###", multicell.Norm2(gaxis), multicell.Norm2(paxis))
		log.Printf("#\t%e\t%e\t%e\t%e", x0, y0, x1, y1)
	}
	log.Printf("Dumping start")
	for gen := 1; gen <= epochlength; gen++ {
		ofilename := fmt.Sprintf("%s_%3.3d.dat", PG_Filename, gen)
		fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintf(fout, "#\t Geno+e0     \tPheno0     \tGeno+e1     \tPheno1   ")
		fmt.Fprintf(fout, "\t||p0-e0||  \t||p1-e1||  \tFit     \tWagFit\n")

		jfilename := fmt.Sprintf("%s_%3.3d.json", json_in, gen)
		pop := multicell.NewPopulation(*ncellsP, *maxpopP)
		pop.FromJSON(jfilename)
		gt0 := pop.GetFlatGenome(multicell.IAncEnv)
		gt1 := pop.GetFlatGenome(multicell.INovEnv)
		et0 := pop.GetFlatStateVec("E", 0)
		et1 := pop.GetFlatStateVec("E", 1)
		for k, g := range gt0 {
			et0[k] = append(et0[k], g...)
			et1[k] = append(et1[k], gt1[k]...)
		}

		pt0 := pop.GetFlatStateVec("P", 0)
		pt1 := pop.GetFlatStateVec("P", 1)

		for k := range et0 {
			multicell.DiffVecs(et0[k], et0[k], mg)
			multicell.DiffVecs(et1[k], et1[k], mg)
			multicell.DiffVecs(pt0[k], pt0[k], mp)
			multicell.DiffVecs(pt1[k], pt1[k], mp)

			x0 := multicell.DotVecs(et0[k], gaxis) / multicell.Norm2(et0[k])
			y0 := multicell.DotVecs(pt0[k], paxis) / multicell.Norm2(pt0[k])
			x1 := multicell.DotVecs(et1[k], gaxis) / multicell.Norm2(et1[k])
			y1 := multicell.DotVecs(pt1[k], paxis) / multicell.Norm2(pt1[k])
			fmt.Fprintf(fout, "\t%e\t%e\t%e\t%e", x0, y0, x1, y1)

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
