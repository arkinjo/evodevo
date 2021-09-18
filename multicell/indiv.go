package multicell

import (
	//"fmt"
	"math"
	"math/rand"
)

type Genome struct { //Genome of an individual
	E  Spmat //Environment cue effect on epigenome
	F  Spmat //Regulatory effect of epigenetic markers on gene expression
	G  Spmat //Regulatory effect of gene expression on epigenetic markers
	Hg Spmat //Contribution of gene expression on higher order complexes
	Hh Spmat //Interaction between higher order complexes
	P  Spmat //Resulting expressed phenotype
	Z  Spmat //Gene expression of offspring
}

type Cell struct { //A 'cell' is characterized by its gene expression and phenotype
	E Vec // Environment encountered by cell; id already in cue
	F Vec // Epigenetic markers
	G Vec // Gene expression
	H Vec // Higher order complexes
	P Vec // Phenotype; id already in cue
}

type Cells struct { //Do we want to reimplement this?
	Ctypes []Cell // Array of cells of different types
}

type Indiv struct { //An individual as an unicellular organism
	Id      int
	DadId   int
	MomId   int
	Genome  Genome
	Copies  []Cells //0: No env; 1: Previous env; 2: Current env
	Z       Vec     // Initial gene expression of offspring
	F0      float64 //Fitness without cues
	Fit     float64 //Fitness with cues
	Util    float64 //Fitness Utility of cues
	CuePlas float64 //Cue Plasticity
	ObsPlas float64 //Observed Plasticity
	Pp      float64 //Degree of polyphenism
}

func NewGenome() Genome { //Generate new genome matrix ensemble
	E := NewSpmat(ngenes, nenv+ncells)
	F := NewSpmat(ngenes, ngenes)
	G := NewSpmat(ngenes, ngenes)
	Hg := NewSpmat(ngenes, ngenes)
	Hh := NewSpmat(ngenes, ngenes)
	P := NewSpmat(ngenes, nenv+ncells)
	Z := NewSpmat(ngenes, ngenes)

	genome := Genome{E, F, G, Hg, Hh, P, Z}

	return genome
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
	for _, r := range G.Hg.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.Hh.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.P.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
	for _, r := range G.Z.Mat {
		for j := range r { //range over keys
			delete(r, j)
		}
	}
}

func (parent *Genome) Copy() Genome { //creates copy of genome
	e := parent.E.Copy()
	f := parent.F.Copy()
	g := parent.G.Copy()
	hg := parent.Hg.Copy()
	hh := parent.Hh.Copy()
	p := parent.P.Copy()
	z := parent.Z.Copy()

	genome := Genome{e, f, g, hg, hh, p, z}

	return genome
}

func DiffGenomes(Gout, G1, G0 Genome) { //Elementwise difference between two genomes
	if withCue {
		Gout.E = DiffSpmat(&G1.E, &G0.E)
	}
	if epig {
		Gout.F = DiffSpmat(&G1.F, &G0.F)
	}

	Gout.G = DiffSpmat(&G1.G, &G0.G)

	if hoc {
		Gout.Hg = DiffSpmat(&G1.Hg, &G0.Hg)
		if hoi {
			Gout.Hh = DiffSpmat(&G1.Hh, &G0.Hh)
		}
	}

	Gout.P = DiffSpmat(&G1.P, &G0.P)
	Gout.Z = DiffSpmat(&G1.Z, &G0.Z)
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
		for _, m := range G.Hg.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}

		if hoi {
			for _, m := range G.Hh.Mat {
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
	for _, m := range G.Z.Mat {
		for _, v := range m {
			lambda2 += v * v
		}
	}

	lambda := math.Sqrt(lambda2)

	if withCue {
		for i, m := range eG.E.Mat {
			for j := range m {
				eG.E.Mat[i][j] = eG.E.Mat[i][j] / lambda
			}
		}
	}

	if epig {
		for i, m := range eG.F.Mat {
			for j := range m {
				eG.F.Mat[i][j] = eG.F.Mat[i][j] / lambda
			}
		}
	}

	for i, m := range eG.G.Mat {
		for j := range m {
			eG.G.Mat[i][j] = eG.G.Mat[i][j] / lambda
		}
	}

	if hoc {
		for i, m := range eG.Hg.Mat {
			for j := range m {
				eG.Hg.Mat[i][j] = eG.Hg.Mat[i][j] / lambda
			}
		}
		if hoi {
			for i, m := range eG.Hh.Mat {
				for j := range m {
					eG.Hh.Mat[i][j] = eG.Hh.Mat[i][j] / lambda
				}
			}
		}
	}

	for i, m := range eG.P.Mat {
		for j := range m {
			eG.P.Mat[i][j] = eG.P.Mat[i][j] / lambda
		}
	}
	for i, m := range eG.Z.Mat {
		for j := range m {
			eG.Z.Mat[i][j] = eG.Z.Mat[i][j] / lambda
		}
	}
	return eG
}

func NewCell(id int) Cell { //Creates a new cell given id of cell.
	e := NewCue(nenv, id)
	f := NewVec(ngenes)
	g := NewVec(ngenes)
	h := NewVec(ngenes)
	p := NewCue(nenv, id)

	cell := Cell{e, f, g, h, p}

	return cell
}

func (cell *Cell) Copy() Cell {
	id := GetId(cell.P) //Extract id part of cell phenotype
	cell1 := NewCell(id)
	cell1.E = CopyVec(cell.E)
	cell1.F = CopyVec(cell.F)
	cell1.G = CopyVec(cell.G)
	cell1.H = CopyVec(cell.H)
	cell1.P = CopyVec(cell.P)

	return cell1
}

func NewCells(ncells int) Cells { // Creates an array of new cells of length Ncells
	Clist := make([]Cell, ncells)
	for id := range Clist {
		Clist[id] = NewCell(id) //Initialize each cell
	}
	Cellarray := Cells{Clist}
	return Cellarray
}

func (cells *Cells) Copy() Cells {
	//Ncells = len(cells.Ctypes)
	Clist := cells.Ctypes
	Cells1 := NewCells(ncells)
	Clist1 := Cells1.Ctypes
	for i, cell := range Clist {
		Clist1[i] = cell.Copy()
	}
	return Cells1
}

func NewIndiv(id int) Indiv { //Creates a new individual
	genome := NewGenome()
	cellcopies := make([]Cells, 3)
	for i := range cellcopies {
		cellcopies[i] = NewCells(ncells)
	}

	z := NewVec(ngenes)

	indiv := Indiv{id, 0, 0, genome, cellcopies, z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	return indiv
}

func (indiv *Indiv) Copy() Indiv { //Deep copier
	indiv1 := NewIndiv(indiv.Id)
	indiv1.DadId = indiv.DadId
	indiv1.MomId = indiv.MomId
	indiv1.Genome = indiv.Genome.Copy()
	for i, ccopy := range indiv.Copies {
		indiv1.Copies[i] = ccopy.Copy()
	}
	indiv1.Z = CopyVec(indiv.Z)
	indiv1.F0 = indiv.F0
	indiv1.Fit = indiv.Fit
	indiv1.Util = indiv.Util
	indiv1.CuePlas = indiv.CuePlas
	indiv1.ObsPlas = indiv.ObsPlas
	indiv1.Pp = indiv.Pp

	return indiv1
}

/*
func (indiv *Indiv) ZeroGenome() { //Empties genome of individual
	indiv.Genome.Zero()
}
*/
func (indiv *Indiv) Mutate() { //Mutates portion of genome of an individual
	r := rand.Intn(genelength) //Randomly choose one of the genome matrices to mutate; with prob proportional to no. of columns
	t := 0

	//This version is specialized for current definition of full model.

	if withCue {
		t += nenv + ncells
		if r < t {
			indiv.Genome.E.mutateSpmat()
		}
	}
	if epig {
		t += ngenes
		if r < t {
			indiv.Genome.F.mutateSpmat()
		}
	}
	t += ngenes
	if r < t {
		indiv.Genome.G.mutateSpmat()
	}
	if hoc {
		t += ngenes
		if r < t {
			indiv.Genome.Hg.mutateSpmat()
		}
		if hoi {
			t += ngenes
			if r < t {
				indiv.Genome.Hh.mutateSpmat()
			}
		}
	}
	t += nenv + ncells
	if r < t {
		indiv.Genome.P.mutateSpmat()
	}
	t += ngenes
	if r < t {
		indiv.Genome.Z.mutateSpmat()
	}

	return
}

func Crossover(dadg, momg *Genome, dadz, momz Vec) (Genome, Genome, Vec, Vec) { //Crossover
	ng0 := dadg.Copy()
	ng1 := momg.Copy()
	nz0 := dadz
	nz1 := momz

	for i := 0; i < ngenes; i++ {
		r := rand.Float64()
		if r < 0.5 {
			e := ng0.E.Mat[i]
			f := ng0.F.Mat[i]
			g := ng0.G.Mat[i]
			hg := ng0.Hg.Mat[i]
			hh := ng0.Hh.Mat[i]
			p := ng0.P.Mat[i]
			z := ng0.Z.Mat[i]
			z0 := nz0[i]

			ng0.E.Mat[i] = ng1.E.Mat[i]
			ng0.F.Mat[i] = ng1.F.Mat[i]
			ng0.G.Mat[i] = ng1.G.Mat[i]
			ng0.Hg.Mat[i] = ng1.Hg.Mat[i]
			ng0.Hh.Mat[i] = ng1.Hh.Mat[i]
			ng0.P.Mat[i] = ng1.P.Mat[i]
			ng0.Z.Mat[i] = ng1.Z.Mat[i]
			nz0[i] = nz1[i]

			ng1.E.Mat[i] = e
			ng1.F.Mat[i] = f
			ng1.G.Mat[i] = g
			ng1.Hg.Mat[i] = hg
			ng1.Hh.Mat[i] = hh
			ng1.P.Mat[i] = p
			ng1.Z.Mat[i] = z
			nz1[i] = z0
		}
	}

	return ng0, ng1, nz0, nz1
}

func Mate(dad, mom *Indiv) (Indiv, Indiv) { //Generates offspring
	genome0, genome1, g0, g1 :=
		Crossover(&dad.Genome, &mom.Genome, dad.Z, mom.Z)

	cells0 := make([]Cells, 3)
	for i := range cells0 {
		cells0[i] = NewCells(ncells)
	}
	cells1 := make([]Cells, 3)
	for i := range cells1 {
		cells1[i] = NewCells(ncells)
	}

	kid0 := Indiv{dad.Id, dad.Id, mom.Id, genome0, cells0, g0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	kid1 := Indiv{mom.Id, dad.Id, mom.Id, genome1, cells1, g1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	kid0.Mutate()
	kid1.Mutate()
	//Mutation happens with 100% probability

	return kid0, kid1
}

func (cells *Cells) get_fitness(envs Cues) float64 {
	d2 := 0.0

	//declaring arrays to be filled
	env := make([]float64, nenv+ncells)
	p := make([]float64, nenv+ncells)
	//id := make([]float64, ncells)
	//idp := make([]float64, ncells)
	for i, cell := range cells.Ctypes {
		env = envs[i]
		p = cell.P
		d2 += dist2Vecs(p, env)
	}
	return math.Exp(-selStrength * d2)
}

func (indiv *Indiv) get_cue_plasticity() float64 { //cue plasticity of individual
	d2 := 0.0
	p0 := NewVec(nenv)
	p := NewVec(nenv)
	Copies := indiv.Copies
	for i, cell := range Copies[2].Ctypes {
		p = cell.P
		p0 = Copies[0].Ctypes[i].P
		d2 += dist2Vecs(p, p0)
	}
	d2 = dist2Vecs(p, p0) / float64(nenv*ncells) //Divide by number of phenotypes to normalize
	return d2
}

func (indiv *Indiv) get_obs_plasticity() float64 { //cue plasticity of individual
	d2 := 0.0
	p0 := NewVec(nenv)
	p := NewVec(nenv)
	Copies := indiv.Copies
	for i, cell := range Copies[2].Ctypes {
		p = cell.P
		p0 = Copies[1].Ctypes[i].P
		d2 += dist2Vecs(p, p0)
	}
	d2 = dist2Vecs(p, p0) / float64(nenv*ncells) //Divide by number of phenotypes to normalize
	return d2
}

func (indiv *Indiv) get_vp() float64 { //Get sum of elementwise variance of phenotype
	pvec := make([]Cue, 0)
	for _, c := range indiv.Copies[2].Ctypes {
		pvec = append(pvec, c.P) //Note: To be used AFTER development
	}
	sigma2p := GetCueVar(pvec)
	return sigma2p
}

func (indiv *Indiv) get_pp(envs Cues) float64 { //Degree of polyphenism of individual; normalize with variance of environment cue
	sigma2env := GetCueVar(envs)
	if sigma2env == 0 {
		return 0
	} else {
		sigma2p := indiv.get_vp()
		return sigma2p / sigma2env
	}
}

func (cell *Cell) DevCell(G Genome, g0 Vec, env Cue) Cell { //Develops a cell given cue
	var diff float64

	h0 := make([]float64, ngenes) //No higher order complexes in embryonic stage
	ve := make([]float64, ngenes)
	vf := make([]float64, ngenes)
	vg := g0
	vh := make([]float64, ngenes)
	vp := make([]float64, nenv+ncells)
	f1 := make([]float64, ngenes)
	g1 := make([]float64, ngenes)
	h1 := make([]float64, ngenes)

	for nstep := 0; nstep < MaxDevStep; nstep++ {
		multMatVec(ve, G.E, env)
		multMatVec(vf, G.G, g0)
		if withCue { //Model with or without cues
			addVecs(f1, vf, ve)
		} else {
			f1 = vf
		}
		applyFnVec(sigmaf, f1)
		if epig { //Allow or disallow epigenetic layer
			multMatVec(g1, G.F, f1)
			applyFnVec(sigmag, g1)
		} else { //Remove epigenetic layer if false
			g1 = f1
		}
		if hoc { //If layer for higher order complexes is present
			multMatVec(vg, G.Hg, g1)
			multMatVec(vh, G.Hh, h0)
			if hoi { //If interactions between higher order complexes is present
				addVecs(h1, vg, vh)
			} else {
				h1 = vg
			}
		} else {
			h1 = g1
		}
		applyFnVec(sigmah, h1)
		multMatVec_T(vp, G.P, h1)
		diff = dist2Vecs(vg, g0)
		g0 = g1
		h0 = h1
		if diff < epsDev {
			break
		}

	}
	//fmt.Println("Phenotype after development:",vpc)
	//fmt.Println("Id after development:",vpid)

	cell.E = env
	cell.F = f1
	cell.G = g1
	cell.H = h1
	cell.P = vp

	return *cell
}

func (cells *Cells) DevCells(G Genome, g0 Vec, envs Cues) Cells {
	env := NewCue(nenv, 0)

	for i := range cells.Ctypes {
		env = envs[i]
		cells.Ctypes[i].DevCell(G, g0, env)
	}

	return *cells
}

func (indiv *Indiv) CompareDev(env, env0 Cues) Indiv { //Compare developmental process under different conditions

	devenv := AddNoisetoCues(env, DevNoise)
	devenv0 := AddNoisetoCues(env0, DevNoise)
	selenv := AddNoisetoCues(env, EnvNoise)
	Clist := indiv.Copies
	//ncells := len(Clist[2].Ctypes)
	//nenv := len(Clist[2].Ctypes[0].E)
	zero := NewCues(ncells, nenv+ncells)

	Clist[0].DevCells(indiv.Genome, indiv.Z, zero)    //Develop without cues
	Clist[1].DevCells(indiv.Genome, indiv.Z, devenv0) //Develop in ancestral (previous) environment
	Clist[2].DevCells(indiv.Genome, indiv.Z, devenv)  //Develop in novel (present) environment

	multMatVec_T(indiv.Z, indiv.Genome.Z, Clist[2].Ctypes[0].G)
	indiv.F0 = Clist[0].get_fitness(selenv)  //Fitness without cues
	indiv.Fit = Clist[2].get_fitness(selenv) //Fitness with cues
	indiv.Util = indiv.Fit - indiv.F0

	indiv.CuePlas = indiv.get_cue_plasticity()
	indiv.ObsPlas = indiv.get_obs_plasticity()
	indiv.Pp = indiv.get_pp(devenv)

	return *indiv
}
