package multicell

import (
	//	"errors"
	//	"fmt"
	//	"log"
	"math"
	"math/rand"
)

type Genome struct { //Genome of an individual
	E Spmat //Environment cue effect on epigenome
	F Spmat //Regulatory effect of epigenetic markers on gene expression
	G Spmat //Regulatory effect of gene expression on epigenetic markers
	H Spmat //Contribution of gene expression on higher order complexes
	J Spmat //Interaction between higher order complexes
	P Spmat //Resulting expressed phenotype
	//Z  Spmat //Gene expression of offspring
}

func NewGenome() Genome { //Generate new genome matrix ensemble
	E := NewSpmat(ngenes, nenv+ncells)
	F := NewSpmat(ngenes, ngenes)
	G := NewSpmat(ngenes, ngenes)
	H := NewSpmat(ngenes, ngenes)
	J := NewSpmat(ngenes, ngenes)
	P := NewSpmat(ngenes, nenv+ncells)
	//Z := NewSpmat(ngenes, ngenes)

	genome := Genome{E, F, G, H, J, P}

	return genome
}

func (G *Genome) Randomize() {

	if withCue {
		G.E.Randomize(CueResponseDensity, sdE)
	}

	if epig {
		G.F.Randomize(GenomeDensity, sdF)
	}
	G.G.Randomize(GenomeDensity, sdG)

	if hoc {
		G.H.Randomize(GenomeDensity, sdH)
		if hoi {
			G.J.Randomize(GenomeDensity, sdJ)
		}
	}
	G.P.Randomize(CueResponseDensity, sdP)
}

func (G *Genome) Clear() { //Sets all entries of genome to zero
	for _, r := range G.E.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.F.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.G.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.H.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.J.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.P.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	/*
		for _, r := range G.Z.Mat {
			for j := range r { //range over keys
				delete(r, j)
			}
		}
	*/
}

func (parent *Genome) Copy() Genome { //creates copy of genome
	e := parent.E.Copy()
	f := parent.F.Copy()
	g := parent.G.Copy()
	hg := parent.H.Copy()
	hh := parent.J.Copy()
	p := parent.P.Copy()

	genome := Genome{e, f, g, hg, hh, p}

	return genome
}

func DiffGenomes(Gout, G1, G0 *Genome) { //Elementwise difference between two genomes
	if withCue {
		Gout.E = DiffSpmat(&G1.E, &G0.E)
	}
	if epig {
		Gout.F = DiffSpmat(&G1.F, &G0.F)
	}

	Gout.G = DiffSpmat(&G1.G, &G0.G)

	if hoc {
		Gout.H = DiffSpmat(&G1.H, &G0.H)
		if hoi {
			Gout.J = DiffSpmat(&G1.J, &G0.J)
		}
	}

	Gout.P = DiffSpmat(&G1.P, &G0.P)
}

func (G *Genome) NormalizeGenome() Genome {
	lambda2 := 0.0
	eG := G.Copy()

	if withCue {
		for _, m := range G.E.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}
	}

	if epig {
		for _, m := range G.F.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}
	}

	for _, m := range G.G.Mat {
		for _, v := range m {
			lambda2 += v * v
		}
	}

	if hoc {
		for _, m := range G.H.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}

		if hoi {
			for _, m := range G.J.Mat {
				for _, v := range m {
					lambda2 += v * v
				}
			}
		}
	}

	for _, m := range G.P.Mat {
		for _, v := range m {
			lambda2 += v * v
		}
	}

	if lambda2 == 0 {
		return eG //avoid division by zero
	}

	lambda := math.Sqrt(lambda2)
	sca := 1.0 / lambda
	if withCue {
		eG.E.Scale(sca)
	}

	if epig {
		eG.F.Scale(sca)
	}

	eG.G.Scale(sca)
	if hoc {
		eG.H.Scale(sca)
		if hoi {
			eG.J.Scale(sca)
		}
	}

	eG.P.Scale(sca)

	return eG
}

func Crossover(dadg, momg *Genome) (Genome, Genome) { //Crossover
	ng0 := dadg.Copy()
	ng1 := momg.Copy()
	//nz0 := dadz
	//nz1 := momz

	for i := 0; i < ngenes; i++ {
		r := rand.Float64()
		if r < 0.5 {
			e := ng0.E.Mat[i]
			f := ng0.F.Mat[i]
			g := ng0.G.Mat[i]
			hg := ng0.H.Mat[i]
			hh := ng0.J.Mat[i]
			p := ng0.P.Mat[i]

			ng0.E.Mat[i] = ng1.E.Mat[i]
			ng0.F.Mat[i] = ng1.F.Mat[i]
			ng0.G.Mat[i] = ng1.G.Mat[i]
			ng0.H.Mat[i] = ng1.H.Mat[i]
			ng0.J.Mat[i] = ng1.J.Mat[i]
			ng0.P.Mat[i] = ng1.P.Mat[i]

			ng1.E.Mat[i] = e
			ng1.F.Mat[i] = f
			ng1.G.Mat[i] = g
			ng1.H.Mat[i] = hg
			ng1.J.Mat[i] = hh
			ng1.P.Mat[i] = p
		}
	}

	return ng0, ng1
}
