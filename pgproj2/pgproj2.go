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
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	refgenPtr := flag.Int("refgen", 50, "reference generation for evolved genotype")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "JSON file of input population") //default to empty string

	flag.Parse()

	epochlength := *genPtr
	fmt.Println("epochlength", epochlength)
	refgen := *refgenPtr
	PG_Filename = *pgfilenamePtr

	json_in = *jsoninPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename := fmt.Sprintf("%s_001.json", json_in)
	pop0.FromJSON(jfilename)
	multicell.SetParams(pop0.Params)

	pop1 := multicell.NewPopulation(*ncellsP, *maxpopP) //with randomized genome to start
	jfilename = fmt.Sprintf("%s_%3.3d.json", json_in, refgen)
	fmt.Println("Reference population :", jfilename)
	pop1.FromJSON(jfilename)

	fmt.Println("Initialization of population complete")
	dtint := time.Since(t0)
	fmt.Println("Time taken for initialization : ", dtint)

	fmt.Println("Dumping projections")
	tdump := time.Now()

	// Reference direction
	env0 := multicell.FlattenEnvs(pop0.AncEnvs)
	env1 := multicell.FlattenEnvs(pop0.NovEnvs)
	lenE := len(env0)
	denv := multicell.NewVec(lenE)
	multicell.DiffVecs(denv, env1, env0)
	multicell.NormalizeVec(denv)

	// Get Principal Axis of <phenotype-genotype> cross-covariance
	log.Println("Finding Principal Axes")
	gmix := pop0.GetFlatGenome()
	pmix := pop0.GetFlatStateVec("P", 0)

	g1 := pop1.GetFlatGenome()
	p1 := pop1.GetFlatStateVec("P", 1)

	gmix = append(gmix, g1...)
	pmix = append(pmix, p1...)

	mp, mg, cov := multicell.GetCrossCov(pmix, gmix, true, true)
	U, _, V := multicell.GetSVD(cov)
	paxis := make([]float64, len(mp))
	gaxis := make([]float64, len(mg))
	for i := range paxis {
		paxis[i] = U.At(i, 0)
	}
	for i := range gaxis {
		gaxis[i] = V.At(i, 0)
	}
	neg := multicell.DotVecs(paxis, denv)
	if neg < 0.0 {
		multicell.ScaleVec(paxis, -1.0, paxis)
		multicell.ScaleVec(gaxis, -1.0, gaxis)
	}

	log.Printf("Dumping start")
	for gen := 1; gen <= epochlength; gen++ { //Also project population after pulling back to ancestral environment.
		ofilename := fmt.Sprintf("%s_%3.3d.dat", PG_Filename, gen)
		fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(fout, "#AncPhen \t NovPhen \t Genotype")
		jfilename := fmt.Sprintf("%s_%3.3d.json", json_in, gen)
		pop := multicell.NewPopulation(*ncellsP, *maxpopP)
		pop.FromJSON(jfilename)
		pt0 := pop.GetFlatStateVec("P", 0)
		pt1 := pop.GetFlatStateVec("P", 1)
		gt := pop.GetFlatGenome()
		for k, g := range gt {
			multicell.DiffVecs(gt[k], g, mg)
			multicell.DiffVecs(pt0[k], pt0[k], mp)
			multicell.DiffVecs(pt1[k], pt1[k], mp)

			x0 := multicell.DotVecs(pt0[k], paxis)
			x1 := multicell.DotVecs(pt1[k], paxis)
			y := multicell.DotVecs(gt[k], gaxis)
			fmt.Fprintf(fout, "%f\t %f\t %f\n", x0, x1, y)
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
