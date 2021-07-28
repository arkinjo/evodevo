package unicell

import (
	//"fmt"
	"math"
	"math/rand"
)

type Genome struct { //Genome of an individual
	G Spmat
	E Spmat
	Z Spmat
	P Spmat
}

type Cell struct { //A 'cell' is characterized by its gene expression and phenotype
	G Vec // Gene expression
	P Cue // Phenotype
	E Cue // Environment encountered by cell
}

type Indiv struct { //An individual as an unicellular organism
	Id       int
	DadId	 int
	MomId	 int
	F   	 float64 //Fitness in novel (current) environment
	F0 		 float64 //Fitness in previous environment
	Genome   Genome
	Cells    []Cell
	Z        Vec // Initial gene expression of offspring.
	Pl		 float64 // Degree of observed plasticity
	Plc		 float64 // Degree of cue plasticity
	u		 float64 // Utility of plasticity
}

func NewGenome() Genome { //Generate new genome
	G := NewSpmat(Ngenes, Ngenes, GenomeDensity)
	E := NewSpmat(Ngenes, Nenv, GenomeDensity)
	Z := NewSpmat(Ngenes, Ngenes, GenomeDensity)
	P := NewSpmat(Ngenes, Nenv, GenomeDensity)

	genome := Genome{G, E, Z, P}

	return genome
}

func (parent *Genome) Copy() Genome { //creates copy of parent's genome
	g := make(Spmat, Ngenes)
	e := make(Spmat, Ngenes)
	z := make(Spmat, Ngenes)
	p := make(Spmat, Ngenes)

	for i, m := range parent.G {
		g[i] = make(map[int]float64)
		for j, v := range m {
			g[i][j] = v
		}
	}
	for i, m := range parent.E {
		e[i] = make(map[int]float64)
		for j, v := range m {
			e[i][j] = v
		}
	}
	for i, m := range parent.Z {
		z[i] = make(map[int]float64)
		for j, v := range m {
			z[i][j] = v
		}
	}
	for i, m := range parent.P {
		p[i] = make(map[int]float64)
		for j, v := range m {
			p[i][j] = v
		}
	}

	genome := Genome{g, e, z, p}

	return genome
}

func DiffGenomes(Gout, G0, G1 Genome){ //Elementwise difference between two genomes
	for i := 0 ; i < Ngenes ; i++ {
		for j := 0; j < Ngenes; j++ {
			Gout.G[i][j] = G0.G[i][j] - G1.G[i][j]
		}
		for j := 0; j < Nenv; j++ {
			Gout.E[i][j] = G0.E[i][j] - G1.E[i][j]
		}
		for j := 0; j < Nenv; j++ {
			Gout.P[i][j] = G0.P[i][j] - G1.P[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.Z[i][j] = G0.Z[i][j] - G1.Z[i][j]
		}
	}
}

func(G *Genome) NormalizeGenome() Genome{
	lambda2 := 0.0
	eG := G.Copy()
	
	for _,m := range G.G{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.G{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.G{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.G{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	lambda := math.Sqrt(lambda2)
	for i,m := range eG.G{
		for j := range m{
			eG.G[i][j] = eG.G[i][j]/lambda
		}
	}
	for i,m := range eG.E{
		for j := range m{
			eG.E[i][j] = eG.E[i][j]/lambda
		}
	}
	for i,m := range eG.P{
		for j := range m{
			eG.P[i][j] = eG.P[i][j]/lambda
		}
	}
	for i,m := range eG.Z{
		for j := range m{
			eG.Z[i][j] = eG.Z[i][j]/lambda
		}
	}
	return eG
}

func NewCell() Cell { //Creates a new cell
	g := NewVec(Ngenes)
	p := NewCue(Nenv)
	e := NewCue(Nenv)

	cell := Cell{g, p, e}

	return cell
}

func NewCells(Ncells int) []Cell { // Creates an array of new cells of length Ncells
	Clist := make([]Cell, Ncells)
	for i := range Clist {
		Clist[i] = NewCell() //Initialize each cell
	}
	return Clist
}

func NewIndiv(id int) Indiv { //Creates a new individual
	genome := NewGenome()
	cells := NewCells(Ncells)
	z := NewVec(Ngenes)

	indiv := Indiv{id, 0, 0, 0.0, 0.0, genome, cells, z, 0.0, 0.0, 0.0}

	return indiv
}

func (indiv *Indiv) Mutate() { //Mutates genome of an individual
	ipos := rand.Intn(GeneLength)
	if ipos < Ngenes {
		mutateSpmat(indiv.Genome.G, Ngenes)
	} else if ipos < Ngenes+Nenv {
		mutateSpmat(indiv.Genome.E, Nenv)
	} else if ipos < 2*Ngenes+Nenv {
		mutateSpmat(indiv.Genome.Z, Ngenes)
	} else {
		mutateSpmat(indiv.Genome.P, Nenv)
	}

	return
}

func Crossover(dadg, momg *Genome, dadz, momz Vec) (Genome, Genome, Vec, Vec) { //Crossover
	ng0 := dadg.Copy()
	ng1 := momg.Copy()
	nz0 := dadz
	nz1 := momz

	for i := 0; i < Ngenes; i++ {
		r := rand.Float64()
		if r < 0.5 {
			g := ng0.G[i]
			e := ng0.E[i]
			z := ng0.Z[i]
			p := ng0.P[i]
			z0 := nz0[i]

			ng0.G[i] = ng1.G[i]
			ng0.E[i] = ng1.E[i]
			ng0.Z[i] = ng1.Z[i]
			ng0.P[i] = ng1.P[i]
			nz0[i] = nz1[i]

			ng1.G[i] = g
			ng1.E[i] = e
			ng1.Z[i] = z
			ng1.P[i] = p
			nz1[i] = z0
		}
	}

	return ng0, ng1, nz0, nz1
}

func Mate(dad, mom *Indiv) (Indiv, Indiv) { //Generates offspring
	genome0, genome1, g0, g1 :=
		Crossover(&dad.Genome, &mom.Genome, dad.Z, mom.Z)

	cells0 := NewCells(Ncells)
	cells1 := NewCells(Ncells)

	kid0 := Indiv{dad.Id, dad.Id, mom.Id, 0.0, 0.0, genome0,  cells0, g0, 0.0, 0.0, 0.0}
	kid1 := Indiv{mom.Id, dad.Id, mom.Id, 0.0, 0.0, genome1,  cells1, g1, 0.0, 0.0, 0.0}

	kid0.Mutate()
	kid1.Mutate()

	return kid0, kid1
}

func get_fitness(p, env Vec) float64 { //Difference between phenotype and selection environment
	d2 := dist2Vecs(p, env) //Respect Hamming dist = Euclidean dist squared for Boolean vectors
	return math.Exp(-s * d2)
}

func (indiv *Indiv) get_cue_plasticity() float64 { //cue plasticity of individual
	Clist := indiv.Cells
	p0 := Clist[0].P.C
	p := Clist[2].P.C
	d2 := dist2Vecs(p,p0)/float64(Nenv) //Divide by number of phenotypes to normalize
	return d2
}

func(indiv *Indiv) get_obs_plasticity() float64 { //observed plasticity of individual
	Clist := indiv.Cells
	p0 := Clist[1].P.C
	p := Clist[2].P.C
	d2 := dist2Vecs(p,p0)/float64(Nenv) //Divide by number of phenotypes to normalize
	return d2
}

func (indiv *Indiv) get_utility() float64 { //Fitness utility of plasticity of individual
	F := indiv.F
	F0 := indiv.F0	
	return F-F0
}

func (indiv *Indiv) Develop(env, env0 Cue) Indiv { //Developmental process
	var diff, fness0, fness float64

	devenv := env.AddNoise(DevNoise)
	selenv := env.AddNoise(EnvNoise)
	Clist := indiv.Cells

	Clist[0].G = indiv.Z
	g0 := Clist[0].G
	p0 := Clist[0].P.C

	vg := make(Vec, Ngenes)
	ve := make(Vec, Ngenes)
	g1 := make(Vec, Ngenes)

	for nstep := 0; nstep < MaxDevStep; nstep++ { // Measure without cue; for comparison
		multMatVec(g1, indiv.Genome.G, g0)
		g1 = vg
		applyFnVec(sigma, g1)

		multMatVec_T(p0, indiv.Genome.P, g1)
		applyFnVec(rho, p0)
		fness0 = get_fitness(p0, selenv.C)
		diff = distVecs(g1, g0)
		g0 = g1
		if diff < epsDev { //Convergence criterion
			break
		}
	}

	indiv.F0 = fness0

	Clist[1].G = indiv.Z
	g0 = Clist[1].G
	p0 = Clist[1].P.C

	Clist[1].G = g0
	Clist[1].P.C = p0

	for nstep := 0; nstep < MaxDevStep; nstep++ { // Measure relative to previous epoch; for comparison
		multMatVec(g1, indiv.Genome.G, g0)
		multMatVec(ve, indiv.Genome.E, env0.C) 
		if WithCue {
			addVecs(g1, vg, ve)
		} else {
			g1 = vg
		}
		applyFnVec(sigma, g1)

		multMatVec_T(p0, indiv.Genome.P, g1)
		applyFnVec(rho, p0)
		//fness0 = get_fitness(p0, selenv.C)
		diff = distVecs(g1, g0)
		g0 = g1
		if diff < epsDev { //Convergence criterion
			break
		}
	}

	Clist[1].G = g0
	Clist[1].P.C = p0

	Clist[2].G = indiv.Z
	g := Clist[2].G
	p := Clist[2].P.C

	for nstep := 0; nstep < MaxDevStep; nstep++ { //development with environment cue
		multMatVec(vg, indiv.Genome.G, g)
		multMatVec(ve, indiv.Genome.E, devenv.C)
		if WithCue {
			addVecs(g1, vg, ve)
		} else {
			g1 = vg
		}

		applyFnVec(sigma, g1)

		multMatVec_T(p, indiv.Genome.P, g1)
		applyFnVec(rho, p)
		fness = get_fitness(p, selenv.C)
		diff := distVecs(g1, g)
		g = g1
		if diff < epsDev { //Convergence criterion
			break
		}
	}

	indiv.F = fness

	indiv.Cells[2].G = g
	indiv.Cells[2].P.C = p

	multMatVec_T(indiv.Z, indiv.Genome.Z, g) //Genome of offspring

	indiv.Pl = indiv.get_obs_plasticity()
	indiv.Plc = indiv.get_cue_plasticity()
	indiv.u = indiv.get_utility()

	return *indiv
}
