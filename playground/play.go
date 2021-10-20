package main

import (
	"fmt"
	"flag"
	//	"fmt"
	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"
var jfilename string

func main() {
	seedPtr := flag.Int("seed", 1, "random seed")
	ncelltypesPtr := flag.Int("celltypes", 1, "number of cell types/phenotypes simultaneously trained") //default to unicellular case
	cuestrengthPtr := flag.Float64("cuestrength", 1.0, "control size of var contribution of environmental cue")
	hoistrengthPtr := flag.Float64("hoistrength", 1.0, "control size of var contribution of higher order interactions")
	epigPtr := flag.Bool("epig", true, "Add layer representing epigenetic markers")
	HOCPtr := flag.Bool("HOC", true, "Add layer representing higher order complexes")
	omegaPtr := flag.Float64("omega", 1.0, "parameter of sigmoid")
	nsamplePtr := flag.Int("nsample", 100, "number of samples")
	flag.Parse()

	multicell.SetSeed(int64(*seedPtr))
	multicell.Omega = *omegaPtr
	multicell.SetNcells(*ncelltypesPtr)
	multicell.SetLayers(*cuestrengthPtr, *hoistrengthPtr, *epigPtr, *HOCPtr)

	nsamples := *nsamplePtr
	
	env0 := multicell.RandomEnvs(1, multicell.GetNenv(), 0.5)
	//	env1 := multicell.RandomEnvs(1, multicell.GetNenv(), 0.5)
	for i := 0; i < nsamples; i++ {
		indiv0 := multicell.NewIndiv(i)
		conv0 := false
		conv1 := false
		indiv0.Genome.Randomize()
		_,err0 := indiv0.Copies[0].DevCells(indiv0.Genome, env0)
		_,err1 := indiv0.Copies[1].DevCells(indiv0.Genome, env0)
		if err0 == nil {
			conv0 = true
		}
		if err1 == nil {
			conv1 = true
		}
		dp := multicell.Dist2Vecs(indiv0.Copies[0].Ctypes[0].P,
			indiv0.Copies[1].Ctypes[0].P)
		//		if dp > 1.0e-2 {
			//			fmt.Println("ID:", i, "P0:",indiv0.Copies[0].Ctypes[0].P)
			//			fmt.Println("ID:", i, "P1:",indiv0.Copies[1].Ctypes[0].P)
		fmt.Println("ID:", i, "Diff:", dp, conv0, conv1) 
		//		}

	}
	
}
