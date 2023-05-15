package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"github.com/arkinjo/evodevo/multicell"
)

var PG_Filename string //Dump for phenotypes and genotypes
var json_in string     //JSON encoding of initial population; default to empty string

func main() {
	//Given population, output novel phenotype values
	log.Println("Starting...")
	t0 := time.Now()
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	maxdevstepP := flag.Int("maxdevstep", 200, "maximum number of steps for development")
	ngenesP := flag.Int("ngenes", 200, "number of genes")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")
	//genPtr := flag.Int("ngen", 200, "number of generation/epoch")
	//envP := flag.Bool("env", true, "with environment?")
	//ref1Ptr := flag.String("ref1", "", "reference JSON file 1")
	//ref2Ptr := flag.String("ref2", "", "reference JSON file 2")

	pgfilenamePtr := flag.String("PG_file", "phenogeno", "Filename of projected phenotypes and genotypes")
	jsongzinPtr := flag.String("jsongzin", "", "basename of JSON files")

	flag.Parse()

	settings := multicell.CurrentSettings()
	settings.MaxPop = *maxpopP
	settings.MaxDevStep = *maxdevstepP
	settings.NGenes = *ngenesP
	settings.NCells = *ncellsP
	//epochlength := *genPtr
	//refgen1 := *ref1Ptr
	//refgen2 := *ref2Ptr
	//withEnv := *envP

	PG_Filename = *pgfilenamePtr

	json_in = *jsongzinPtr

	if json_in == "" {
		log.Fatal("Must specify JSON input file.")
	}
	/*
		if refgen1 == "" {
			log.Fatal("Must specify JSON reference file 1.")
		}

			if refgen2 == "" {
				log.Fatal("Must specify JSON reference file 2.")
			}


			log.Println("epochlength", epochlength)
	*/

	tdump := time.Now()

	/*
		log.Println("Reading Pop0")
		pop0 := multicell.NewPopulation(settings)
		//fmt.Println("Reference population :", refgen1)
		//pop0.ImportPopGz(refgen1)
		settings = pop0.Params
		multicell.SetParams(settings)

	*/
	//Nenv := multicell.GetNenv()
	//Nsel := multicell.GetNsel()

	// Reference direction (Selective Envs only)
	///env0 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop0.AncEnvs))
	//env1 := multicell.FlattenEnvs(multicell.GetSelEnvs(pop0.NovEnvs))
	//lenP := len(env0)

	/*
		g00 := pop0.GetFlatGenome(multicell.IAncEnv)
		e00 := pop0.GetFlatStateVec("E", multicell.IAncEnv, 0, Nenv)
		if withEnv {
			for k, e := range e00 {
				g00[k] = append(g00[k], e...)
			}
		}
		log.Println("Reading Pop1")
		pop1 := multicell.NewPopulation(settings)
		fmt.Println("Reference population 2:", refgen2)
		pop1.ImportPopGz(refgen2)

		g11 := pop1.GetFlatGenome(multicell.INovEnv)
		e11 := pop1.GetFlatStateVec("E", multicell.INovEnv, 0, Nenv)
		if withEnv {
			for k, e := range e11 {
				g11[k] = append(g11[k], e...)
			}
		}

		log.Println("Finding Principal Axes")

		mg0 := multicell.GetMeanVec(g00)
		mg1 := multicell.GetMeanVec(g11)
		lenG := len(mg0)
		midg := mg0
		//	midg := multicell.NewVec(lenG)
		//	multicell.AddVecs(midg, mg0, mg1)
		//	multicell.ScaleVec(midg, 0.5, midg)
		gaxis := multicell.NewVec(lenG)
		multicell.DiffVecs(gaxis, mg1, mg0)
		//multicell.NormalizeVec(gaxis)
		len2_gaxis := multicell.Norm2Sq(gaxis)

		midp := env0
		//	midp := multicell.NewVec(lenP)
		//	multicell.AddVecs(midp, env0, env1)
		//	multicell.ScaleVec(midp, 0.5, midp)
		paxis := multicell.NewVec(lenP)
		multicell.DiffVecs(paxis, env1, env0)
		//	multicell.NormalizeVec(paxis)
		len2_paxis := multicell.Norm2Sq(paxis)
	*/

	fmt.Println("Importing population...")
	jfilename := fmt.Sprintf("%s.json.gz", json_in)
	pop := multicell.NewPopulation(settings)
	pop.ImportPopGz(jfilename)
	settings = pop.Params
	multicell.SetParams(settings)
	//Nsel := multicell.GetNsel()

	fmt.Println("Successfully imported population from %s", jfilename)

	ofilename := fmt.Sprintf("%s.phen", *jsongzinPtr)
	fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatal(err)
	}

	maxWagFit := 0.0 //Track maximum rescaled fitness so far
	IndivCount := 0

	for _, indiv := range pop.Indivs {
		if indiv.WagFit > maxWagFit {
			maxWagFit = indiv.WagFit
			IndivCount += 1
			fmt.Fprintf(fout, "Id: %d \t Err: %e \t WagFit: %e \n", indiv.Id, indiv.Dp1e1, indiv.WagFit)
			fmt.Printf("Id: %d \t Err: %e \t WagFit: %e \n", indiv.Id, indiv.Dp1e1, indiv.WagFit)
			for i, body := range indiv.Copy().Bodies {
				fmt.Fprintf(fout, "Body: %d \n", i)
				for j, cell := range body.Copy().Cells {
					fmt.Fprintf(fout, "Cell: %d \n", j)
					for k, np := range cell.P {
						fmt.Fprintf(fout, "Trait: %d \t Phen: %e \t Cue: %e \n", k, np, cell.E[k])
					}
				}
			}
		}

		if maxWagFit > 0.95 || IndivCount > 10 { //Early quick sampling termination
			fmt.Printf("Sampled %d individuals \n", IndivCount)
			fmt.Println("Maximum samples reached")
			break
		}
	}
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	/*
		log.Printf("Dumping start")
		for gen := 1; gen <= epochlength; gen++ {
			ofilename := fmt.Sprintf("%s_%3.3d.dat", PG_Filename, gen)
			fout, err := os.OpenFile(ofilename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
			if err != nil {
				log.Fatal(err)
			}
			fmt.Fprintf(fout, "#\t Geno+e0     \tPheno0     \tGeno+e1     \tPheno1   ")
			fmt.Fprintf(fout, "\t||p0-e0||  \t||p1-e1||  \tFit     \tWagFit\n")

			jfilename := fmt.Sprintf("%s_%3.3d.json.gz", json_in, gen)
			pop := multicell.NewPopulation(settings)
			pop.ImportPopGz(jfilename)
			gt0 := pop.GetFlatGenome(0)
			gt1 := pop.GetFlatGenome(1)
			et0 := pop.GetFlatStateVec("E", 0, 0, Nenv)
			et1 := pop.GetFlatStateVec("E", 1, 0, Nenv)
			if withEnv {
				for k, e := range et0 {
					gt0[k] = append(gt0[k], e...)
					gt1[k] = append(gt1[k], et1[k]...)
				}
			}
			pt0 := pop.GetFlatStateVec("P", 0, 0, Nsel)
			pt1 := pop.GetFlatStateVec("P", 1, 0, Nsel)
			tx0 := multicell.NewVec(lenG)
			tx1 := multicell.NewVec(lenG)
			ty0 := multicell.NewVec(lenP)
			ty1 := multicell.NewVec(lenP)

			dx0 := make([]float64, 0)
			dy0 := make([]float64, 0)
			dx1 := make([]float64, 0)
			dy1 := make([]float64, 0)
			for k := range et0 {
				multicell.DiffVecs(tx0, gt0[k], midg)
				multicell.DiffVecs(tx1, gt1[k], midg)

				multicell.DiffVecs(ty0, pt0[k], midp)
				multicell.DiffVecs(ty1, pt1[k], midp)

				x0 := multicell.DotVecs(tx0, gaxis) / len2_gaxis
				y0 := multicell.DotVecs(ty0, paxis) / len2_paxis

				x1 := multicell.DotVecs(tx1, gaxis) / len2_gaxis
				y1 := multicell.DotVecs(ty1, paxis) / len2_paxis

				fmt.Fprintf(fout, "Data\t%e\t%e\t%e\t%e", x0, y0, x1, y1)

				dx0 = append(dx0, x0)
				dy0 = append(dy0, y0)
				dx1 = append(dx1, x1)
				dy1 = append(dy1, y1)

				dp1e1 := pop.Indivs[k].Dp1e1
				dp0e0 := pop.Indivs[k].Dp0e0
				fit := pop.Indivs[k].Fit
				wf := pop.Indivs[k].WagFit

				fmt.Fprintf(fout, "\t%e\t%e\t%e\t%e\n", dp0e0, dp1e1, fit, wf)

			}
			ax0, sx0 := avesd(dx0)
			ay0, sy0 := avesd(dy0)
			fmt.Fprintf(fout, "AVESD0\t%d\t%e\t%e\t%e\t%e\n", gen, ax0, ay0, sx0, sy0)
			ax1, sx1 := avesd(dx1)
			ay1, sy1 := avesd(dy1)
			fmt.Fprintf(fout, "AVESD1\t%d\t%e\t%e\t%e\t%e\n", gen, ax1, ay1, sx1, sy1)
			err = fout.Close()
			if err != nil {
				log.Fatal(err)
			}

		}

	*/
	dtdump := time.Since(tdump)
	fmt.Println("Time taken :", dtdump)
	//fmt.Printf("Projections written to %s_*.dat \n", PG_Filename)
	dt := time.Since(t0)

	fmt.Println("Total time taken : ", dt)
}
