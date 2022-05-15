package multicell

import (
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
}

func NewGenome() Genome { //Generate new genome matrix ensemble
	E := NewSpmat(ngenes, nenv+ncells)
	F := NewSpmat(ngenes, ngenes)
	G := NewSpmat(ngenes, ngenes)
	H := NewSpmat(ngenes, ngenes)
	J := NewSpmat(ngenes, ngenes)
	P := NewSpmat(nenv+ncells, ngenes)
	genome := Genome{E, F, G, H, J, P}

	return genome
}

func (G *Genome) Randomize() {
	G.E.Randomize(DensityE)
	G.F.Randomize(DensityF)
	G.G.Randomize(DensityF)
	G.H.Randomize(DensityH)
	G.J.Randomize(DensityJ)
	G.P.Randomize(DensityP)
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
	if withE {
		Gout.E = DiffSpmat(&G1.E, &G0.E)
	}
	if withF {
		Gout.F = DiffSpmat(&G1.F, &G0.F)
	}

	Gout.G = DiffSpmat(&G1.G, &G0.G)

	if withH {
		Gout.H = DiffSpmat(&G1.H, &G0.H)
		if withJ {
			Gout.J = DiffSpmat(&G1.J, &G0.J)
		}
	}

	Gout.P = DiffSpmat(&G1.P, &G0.P)
}

func (G *Genome) NormalizeGenome() Genome {
	lambda2 := 0.0
	eG := G.Copy()

	if withE {
		for _, m := range G.E.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}
	}

	if withF {
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

	if withH {
		for _, m := range G.H.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}

		if withJ {
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
	if withE {
		eG.E.Scale(sca)
	}

	if withF {
		eG.F.Scale(sca)
	}

	eG.G.Scale(sca)
	if withH {
		eG.H.Scale(sca)
		if withJ {
			eG.J.Scale(sca)
		}
	}

	eG.P.Scale(sca)

	return eG
}

func (genome *Genome) FlatVec() Vec {
	vec := make([]float64, 0)

	nenvcell := nenv + ncells
	if withE {
		for _, v := range genome.E.Mat {
			for j := 0; j < nenvcell; j++ {
				vec = append(vec, v[j])
			}
		}
	}

	if withF {
		for _, v := range genome.F.Mat {
			for j := 0; j < ngenes; j++ {
				vec = append(vec, v[j])
			}
		}
	}

	for _, v := range genome.G.Mat {
		for j := 0; j < ngenes; j++ {
			vec = append(vec, v[j])
		}
	}

	if withH {
		for _, v := range genome.H.Mat {
			for j := 0; j < ngenes; j++ {
				vec = append(vec, v[j])
			}
		}
		if withJ {
			for _, v := range genome.J.Mat {
				for j := 0; j < ngenes; j++ {
					vec = append(vec, v[j])
				}
			}
		}
	}
	for _, v := range genome.P.Mat {
		for j := 0; j < nenvcell; j++ {
			vec = append(vec, v[j])
		}
	}

	return vec
}

func (genome *Genome) Mutate() {
	r := rand.Intn(geneLength)
	//Randomly choose one of the genome matrices to mutate;
	//with prob proportional to no. of columns
	t := 0

	//This version is specialized for current definition of full model.

	if withE {
		t += nenv + ncells
		if r < t {
			genome.E.mutateSpmat(DensityE, mutRate)
			return
		}
	}
	if withF {
		t += ngenes
		if r < t {
			genome.F.mutateSpmat(DensityF, mutRate)
			return
		}
	}
	t += ngenes
	if r < t {
		genome.G.mutateSpmat(DensityG, mutRate)
		return
	}
	if withH {
		t += ngenes
		if r < t {
			genome.H.mutateSpmat(DensityH, mutRate)
			return
		}
		if withJ {
			t += ngenes
			if r < t {
				genome.J.mutateSpmat(DensityJ, mutRate)
				return
			}
		}
	}
	t += nenv + ncells
	if r < t {
		genome.P.mutateSpmat(DensityP, mutRate)
		return
	}
}
