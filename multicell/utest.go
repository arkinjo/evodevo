package multicell

import (
	"math"
)

func TestEqualGenome(G0, G1 Genome) bool { //Elementwise difference between two genomes
	var d float64

	for i := 0; i < ngenes; i++ {

		if withCue {
			for j := 0; j < nenv+ncells; j++ {
				d += math.Abs(G1.E[i][j] - G0.E[i][j])
			}
		}

		if epig {
			for j := 0; j < ngenes; j++ {
				d += math.Abs(G1.F[i][j] - G0.F[i][j])
			}
		}

		for j := 0; j < ngenes; j++ {
			d += math.Abs(G1.G[i][j] - G0.G[i][j])
		}

		if hoc {
			for j := 0; j < ngenes; j++ {
				d += math.Abs(G1.Hg[i][j] - G0.Hg[i][j])
			}
			if hoi {
				for j := 0; j < ngenes; j++ {
					d += math.Abs(G1.Hh[i][j] - G0.Hh[i][j])
				}
			}
		}

		for j := 0; j < nenv+ncells; j++ {
			d += math.Abs(G1.P[i][j] - G0.P[i][j])
		}
		for j := 0; j < ngenes; j++ {
			d += math.Abs(G1.Z[i][j] - G0.Z[i][j])
		}
	}
	return (d == 0.0)

}
