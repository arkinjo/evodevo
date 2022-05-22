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
	ngenesP := flag.Int("ngenes", 200, "number of genes")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	ref1Ptr := flag.String("ref1", "", "reference JSON file 1")
	ref2Ptr := flag.String("ref2", "", "reference JSON file 2")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsoninPtr := flag.String("jsonin", "", "basename of JSON files")

	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	settings.NGenes = *ngenesP
	settings.NCells = *ncellsP
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

	log.Println("Reading Pop0")
	pop0 := multicell.NewPopulation(settings)
	fmt.Println("Reference population :", refgen1)
	pop0.FromJSON(refgen1)
	settings = pop0.Params
	multicell.SetParams(settings)

	Nsel := multicell.GetNsel()
	Nenv := multicell.GetNenv()

	// Reference direction
	env0 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop0.AncEnvs))
	env1 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop0.NovEnvs))
	lenP := len(env0)
	denv := multicell.NewVec(lenP)
	multicell.DiffVecs(denv, env1, env0)

	g00 := pop0.GetFlatGenome(multicell.IAncEnv)
	e00 := pop0.GetFlatStateVec("E", 0, 0, Nenv)
	for k, g := range g00 {
		e00[k] = append(e00[k], g...)
	}

	log.Println("Reading Pop1")
	pop1 := multicell.NewPopulation(settings)
	fmt.Println("Reference population 2:", refgen2)
	pop1.FromJSON(refgen2)

	g11 := pop1.GetFlatGenome(multicell.INovEnv)
	e11 := pop1.GetFlatStateVec("E", 1, 0, Nenv)
	for k, g := range g11 {
		e11[k] = append(e11[k], g...)
	}

	log.Println("Finding Principal Axes")
	midp := multicell.NewVec(lenP)
	multicell.AddVecs(midp, env0, env1)
	multicell.ScaleVec(midp, 0.5, midp)
	paxis := multicell.CopyVec(denv)
	multicell.NormalizeVec(paxis)

	mg0 := multicell.GetMeanVec(e00)
	mg1 := multicell.GetMeanVec(e11)
	lenG := len(mg0)
	midg := multicell.NewVec(lenG)
	multicell.AddVecs(midg, mg0, mg1)
	multicell.ScaleVec(midg, 0.5, midg)

	gaxis := multicell.NewVec(lenG)
	multicell.DiffVecs(gaxis, mg1, mg0)
	multicell.NormalizeVec(gaxis)

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
		pop := multicell.NewPopulation(settings)
		pop.FromJSON(jfilename)
		gt0 := pop.GetFlatGenome(0)
		gt1 := pop.GetFlatGenome(1)
		et0 := pop.GetFlatStateVec("E", 0, 0, Nenv)
		et1 := pop.GetFlatStateVec("E", 1, 0, Nenv)
		for k, g := range gt0 {
			et0[k] = append(et0[k], g...)
			et1[k] = append(et1[k], gt1[k]...)
		}

		pt0 := pop.GetFlatStateVec("P", 0, 0, Nsel)
		pt1 := pop.GetFlatStateVec("P", 1, 0, Nsel)

		tx0 := multicell.NewVec(lenG)
		tx1 := multicell.NewVec(lenG)
		ty0 := multicell.NewVec(lenP)
		ty1 := multicell.NewVec(lenP)

		for k := range et0 {
			multicell.DiffVecs(tx0, et0[k], midg)
			multicell.DiffVecs(tx1, et1[k], midg)

			multicell.DiffVecs(ty0, pt0[k], midp)
			multicell.DiffVecs(ty1, pt1[k], midp)

			x0 := multicell.DotVecs(tx0, gaxis) / multicell.Norm2(tx0)
			y0 := multicell.DotVecs(ty0, paxis) / multicell.Norm2(ty0)

			x1 := multicell.DotVecs(tx1, gaxis) / multicell.Norm2(tx1)
			y1 := multicell.DotVecs(ty1, paxis) / multicell.Norm2(ty1)

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
