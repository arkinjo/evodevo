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
	ref1Ptr := flag.Int("ref1", 200, "reference generation for evolved genotype")
	ref2Ptr := flag.Int("ref2", -1, "reference generation for evolved genotype")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "JSON file of input population") //default to empty string

	flag.Parse()

	epochlength := *genPtr
	refgen1 := *ref1Ptr
	refgen2 := *ref2Ptr
	PG_Filename = *pgfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}

	log.Println("epochlength", epochlength)

	log.Println("Reading Pop0")
	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP)
	jfilename := fmt.Sprintf("%s_001.json", json_in)
	pop0.FromJSON(jfilename)
	multicell.SetParams(pop0.Params)

	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	tdump := time.Now()

	// Reference direction
	env0 := multicell.FlattenEnvs(pop0.AncEnvs)
	env1 := multicell.FlattenEnvs(pop0.NovEnvs)
	lenE := len(env0)
	denv := multicell.NewVec(lenE)
	multicell.DiffVecs(denv, env1, env0)
	multicell.NormalizeVec(denv)

	// Get Principal Axis of <phenotype-genotype> cross-covariance
	gmix := make([][]float64, 0)
	pmix := make([][]float64, 0)

	g0 := pop0.GetFlatGenome()
	e00 := pop0.GetFlatStateVec("E", 0)
	//	e01 := pop0.GetFlatStateVec("E", 1)
	for k, g := range g0 {
		e00[k] = append(e00[k], g...)
		//		e01[k] = append(e01[k], g...)
	}
	p00 := pop0.GetFlatStateVec("P", 0)
	//	p01 := pop0.GetFlatStateVec("P", 1)
	gmix = append(gmix, e00...)
	pmix = append(pmix, p00...)
	//	gmix = append(gmix, e01...)
	//	pmix = append(pmix, p01...)
	log.Println("lenG,lenP=", len(gmix[0]), len(pmix[0]))

	log.Println("Reading Pop1")
	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP)
	jfilename = fmt.Sprintf("%s_%3.3d.json", json_in, refgen1)
	fmt.Println("Reference population :", jfilename)
	pop1.FromJSON(jfilename)

	g1 := pop1.GetFlatGenome()
	//	e10 := pop1.GetFlatStateVec("E", 0)
	e11 := pop1.GetFlatStateVec("E", 1)
	for k, g := range g1 {
		//		e10[k] = append(e10[k], g...)
		e11[k] = append(e11[k], g...)
	}
	//	p10 := pop1.GetFlatStateVec("P", 1)
	p11 := pop1.GetFlatStateVec("P", 1)
	//	gmix = append(gmix, e10...)
	//	pmix = append(pmix, p10...)
	gmix = append(gmix, e11...)
	pmix = append(pmix, p11...)

	if refgen2 > 0 {
		log.Println("Reading Pop2")
		pop2 := multicell.NewPopulation(*ncellsP, *maxpopP)
		jfilename = fmt.Sprintf("%s_%3.3d.json", json_in, refgen2)
		fmt.Println("Reference population :", jfilename)
		pop2.FromJSON(jfilename)

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
	for gen := 1; gen <= epochlength; gen++ { //Also project population after pulling back to ancestral environment.
		ofilename := fmt.Sprintf("%s_%3.3d.dat", PG_Filename, gen)
		fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintf(fout, "#Geno+e0(0) \tPheno0(0) \tGeno+e1(0) \tPheno1(0)")
		fmt.Fprintf(fout, "\t Geno+e0(1) \tPheno0(1) \tGeno+e1(1)\tPheno1(1)")
		fmt.Fprintf(fout, "\t ||p0-e0|| \t||p1-e1||\n")

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
			dp1e1 := pop.Indivs[k].Dp1e1
			dp0e0 := pop.Indivs[k].Dp0e0

			for a, pa := range paxis {
				x0 := multicell.DotVecs(et0[k], gaxis[a])
				y0 := multicell.DotVecs(pt0[k], pa)
				x1 := multicell.DotVecs(et1[k], gaxis[a])
				y1 := multicell.DotVecs(pt1[k], pa)
				fmt.Fprintf(fout, "%e\t%e\t%e\t%e", x0, y0, x1, y1)
			}
			fmt.Fprintf(fout, "\t%e\t%e\n", dp1e1, dp0e0)

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
