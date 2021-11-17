package multicell

import (
	"errors"
	"fmt"
	"math"
	//	"sort"
)

func TestEqualGenomes(G0, G1 Genome) float64 { //Elementwise difference between two genomes
	var d float64

	for i := 0; i < ngenes; i++ {

		if withCue {
			for j := 0; j < nenv+ncells; j++ {
				d += math.Abs(G1.E.Mat[i][j] - G0.E.Mat[i][j])
			}
		}

		if epig {
			for j := 0; j < ngenes; j++ {
				d += math.Abs(G1.F.Mat[i][j] - G0.F.Mat[i][j])
			}
		}

		for j := 0; j < ngenes; j++ {
			d += math.Abs(G1.G.Mat[i][j] - G0.G.Mat[i][j])
		}

		if hoc {
			for j := 0; j < ngenes; j++ {
				d += math.Abs(G1.Hg.Mat[i][j] - G0.Hg.Mat[i][j])
			}
			if hoi {
				for j := 0; j < ngenes; j++ {
					d += math.Abs(G1.Hh.Mat[i][j] - G0.Hh.Mat[i][j])
				}
			}
		}

		for j := 0; j < nenv+ncells; j++ {
			d += math.Abs(G1.P.Mat[i][j] - G0.P.Mat[i][j])
		}
		/*
			for j := 0; j < ngenes; j++ {
				d += math.Abs(G1.Z.Mat[i][j] - G0.Z.Mat[i][j])
			}
		*/
	}
	return d

}

func TestEqualPopGenomes(pop0, pop1 Population) float64 { //Test for equal genome across individuals in two populations.
	u := 0.0
	pop0.SortPopIndivs()
	pop1.SortPopIndivs() //Sort before comparison

	for k, indiv := range pop0.Indivs {
		u += TestEqualGenomes(indiv.Genome, pop1.Indivs[k].Genome) //Update whether individual wise genomes are same
	}
	return u //Warning! Ordering of population individuals is important.
}

func DeepVec2Test(v1, v2 [][]float64) float64 {
	u := 0.0
	for i, v := range v1 {
		u += DistVecs1(v, v2[i])
	}
	return u
}

func DeepVec3NovTest(v [][]float64, vv [][][]float64) error { //Check if v is already in vv
	var errtext string
	nov := true
	d := 0.0
	ancindices := make([]int, 0)
	for i, u := range vv {
		d = DeepVec2Test(v, u)
		if d == 0 {
			ancindices = append(ancindices, i)
			nov = false
		}
	}
	if nov {
		return nil
	} else {
		errtext = fmt.Sprint("Input vector same as memory indices: ", ancindices)
		return errors.New(errtext)
	}
}

/*
func DeepVec3Test(v1, v2 [][][]float64) float64 {
	u := 0.0
	for i, v := range v1 {
		u += DeepVec2Test(v, v2[i])
	}
	return u
}
*/
