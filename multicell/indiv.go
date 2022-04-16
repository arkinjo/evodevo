package multicell

import (
	"errors"
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

type Cell struct { //A 'cell' is characterized by its gene expression and phenotype
	E        Vec // Epigenetic markers
	F        Vec // Epigenetic markers
	G        Vec // Gene expression
	H        Vec // Higher order complexes
	P        Vec // Phenotype; id already in cue
	NDevStep int // Developmental path length
}

type Body struct { //Do we want to reimplement this?
	Cells []Cell // Array of cells of different types
}

const (
	INoEnv  = iota //No env
	IAncEnv        // Previous env
	INovEnv        // Current env
)

type Indiv struct { //An individual as an unicellular organism
	Id         int
	DadId      int
	MomId      int
	Genome     Genome
	Bodies     []Body  //INoEnv, IAncEnv, INovEnv (see above const.)
	MSE        float64 //Squared error per cue
	Fit        float64 //Fitness with cues
	WagFit     float64 //Wagner relative fitness
	AncCuePlas float64 //Cue Plasticity in ancestral environment
	NovCuePlas float64 //Cue Plasticity in novel environment
	ObsPlas    float64 //Observed Plasticity
	Pp         float64 //Degree of polyphenism
	NDevStep   int     // Maximum number of developmental steps in cells under NovEnv.
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
		G.E.Randomize(CueResponseDensity, sde)
	}

	if epig {
		G.F.Randomize(GenomeDensity, sdf)
	}
	G.G.Randomize(GenomeDensity, sdg)

	if hoc {
		G.H.Randomize(GenomeDensity, sdhg)
		if hoi {
			G.J.Randomize(GenomeDensity, sdhh)
		}
	}
	G.P.Randomize(CueResponseDensity, sdp)
	//G.Z.Randomize(GenomeDensity)
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
	//z := parent.Z.Copy()

	//genome := Genome{e, f, g, hg, hh, p, z}

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
	//Gout.Z = DiffSpmat(&G1.Z, &G0.Z)
	//Remark: Ensure that Gout is initialized and empty before applying operation
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
	/*
		for _, m := range G.Z.Mat {
			for _, v := range m {
				lambda2 += v * v
			}
		}
	*/

	if lambda2 == 0 {
		return eG //avoid division by zero
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
		for i, m := range eG.H.Mat {
			for j := range m {
				eG.H.Mat[i][j] = eG.H.Mat[i][j] / lambda
			}
		}
		if hoi {
			for i, m := range eG.J.Mat {
				for j := range m {
					eG.J.Mat[i][j] = eG.J.Mat[i][j] / lambda
				}
			}
		}
	}

	for i, m := range eG.P.Mat {
		for j := range m {
			eG.P.Mat[i][j] = eG.P.Mat[i][j] / lambda
		}
	}
	/*
		for i, m := range eG.Z.Mat {
			for j := range m {
				eG.Z.Mat[i][j] = eG.Z.Mat[i][j] / lambda
			}
		}
	*/
	return eG
}

func NewCell(id int) Cell { //Creates a new cell given id of cell.
	e := NewCue(nenv, id)
	f := NewVec(ngenes)
	g := NewVec(ngenes)
	h := NewVec(ngenes)
	p := NewCue(nenv, id)
	cell := Cell{e, f, g, h, p, 0}

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
	cell1.NDevStep = cell.NDevStep

	return cell1
}

func NewBody(ncells int) Body { // Creates an array of new cells of length Ncells
	cells := make([]Cell, ncells)
	for id := range cells {
		cells[id] = NewCell(id) //Initialize each cell
	}
	return Body{cells}
}

func (body *Body) Copy() Body {
	body1 := NewBody(ncells)
	for i, cell := range body.Cells {
		body1.Cells[i] = cell.Copy()
	}
	return body1
}

func NewIndiv(id int) Indiv { //Creates a new individual
	genome := NewGenome()
	bodies := make([]Body, 3)
	for i := range bodies {
		bodies[i] = NewBody(ncells)
	}

	mse := 0.0

	indiv := Indiv{id, 0, 0, genome, bodies, mse, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0}

	return indiv
}

func (indiv *Indiv) Copy() Indiv { //Deep copier
	indiv1 := NewIndiv(indiv.Id)
	indiv1.DadId = indiv.DadId
	indiv1.MomId = indiv.MomId
	indiv1.Genome = indiv.Genome.Copy()
	for i, body := range indiv.Bodies {
		indiv1.Bodies[i] = body.Copy()
	}
	indiv1.Fit = indiv.Fit
	indiv1.ObsPlas = indiv.ObsPlas
	indiv1.Pp = indiv.Pp
	indiv1.NDevStep = indiv.NDevStep

	return indiv1
}

func (indiv *Indiv) Mutate() { //Mutates portion of genome of an individual
	r := rand.Intn(genelength) //Randomly choose one of the genome matrices to mutate; with prob proportional to no. of columns
	t := 0

	//This version is specialized for current definition of full model.

	if withCue {
		t += nenv + ncells
		if r < t {
			indiv.Genome.E.mutateSpmat(CueResponseDensity, sde)
		}
	}
	if epig {
		t += ngenes
		if r < t {
			indiv.Genome.F.mutateSpmat(GenomeDensity, sdf)
		}
	}
	t += ngenes
	if r < t {
		indiv.Genome.G.mutateSpmat(GenomeDensity, sdg)
	}
	if hoc {
		t += ngenes
		if r < t {
			indiv.Genome.H.mutateSpmat(GenomeDensity, sdhg)
		}
		if hoi {
			t += ngenes
			if r < t {
				indiv.Genome.J.mutateSpmat(GenomeDensity, sdhh)
			}
		}
	}
	t += nenv + ncells
	if r < t {
		indiv.Genome.P.mutateSpmat(CueResponseDensity, sdp)
	}
	/*
		t += ngenes
		if r < t {
			indiv.Genome.Z.mutateSpmat(GenomeDensity)
		}
	*/

	return
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
			//z := ng0.Z.Mat[i]
			//z0 := nz0[i]

			ng0.E.Mat[i] = ng1.E.Mat[i]
			ng0.F.Mat[i] = ng1.F.Mat[i]
			ng0.G.Mat[i] = ng1.G.Mat[i]
			ng0.H.Mat[i] = ng1.H.Mat[i]
			ng0.J.Mat[i] = ng1.J.Mat[i]
			ng0.P.Mat[i] = ng1.P.Mat[i]
			//ng0.Z.Mat[i] = ng1.Z.Mat[i]
			//nz0[i] = nz1[i]

			ng1.E.Mat[i] = e
			ng1.F.Mat[i] = f
			ng1.G.Mat[i] = g
			ng1.H.Mat[i] = hg
			ng1.J.Mat[i] = hh
			ng1.P.Mat[i] = p
			//ng1.Z.Mat[i] = z
			//nz1[i] = z0
		}
	}

	return ng0, ng1
}

func Mate(dad, mom *Indiv) (Indiv, Indiv) { //Generates offspring
	genome0, genome1 :=
		Crossover(&dad.Genome, &mom.Genome)
	bodies0 := make([]Body, 3)
	for i := range bodies0 {
		bodies0[i] = NewBody(ncells)
	}
	bodies1 := make([]Body, 3)
	for i := range bodies1 {
		bodies1[i] = NewBody(ncells)
	}

	kid0 := Indiv{dad.Id, dad.Id, mom.Id, genome0, bodies0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0}
	kid1 := Indiv{mom.Id, dad.Id, mom.Id, genome1, bodies1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0}

	kid0.Mutate()
	kid1.Mutate()
	//Mutation happens with 100% probability

	return kid0, kid1
}

func (indiv *Indiv) get_sse(envs Cues) float64 { //use after development; returns sum of squared errors
	sse := 0.0

	for i, cell := range indiv.Bodies[INovEnv].Cells { //range over all cells
		sse += Dist2Vecs(cell.P, envs[i])
	}

	return sse
}

func (indiv *Indiv) get_fitness() float64 { //fitness in novel/present environment
	rawfit := math.Exp(-baseSelStrength * indiv.MSE)
	return rawfit
}

func get_plasticity(body0, body1 Body) float64 { //cue plasticity of individual
	d2 := 0.0
	for i, cell := range body0.Cells {
		d2 += Dist2Vecs(cell.P, body1.Cells[i].P)
	}
	//d2 = d2 / float64(nenv*ncells) //Divide by number of phenotypes to normalize
	return d2 / float64(ncells*(ncells+nenv))
}

func (indiv *Indiv) get_varpheno() float64 { //Get sum of elementwise variance of phenotype
	pvec := make([]Cue, 0)
	for _, c := range indiv.Bodies[INovEnv].Cells {
		pvec = append(pvec, c.P) //Note: To be used AFTER development
	}
	sigma2p := GetCueVar(pvec)
	return sigma2p
}

func (indiv *Indiv) get_polyphenism(envs Cues) float64 { //Degree of polyphenism of individual; normalize with variance of environment cue
	sigma2env := GetCueVar(envs)
	if sigma2env == 0 {
		return 0
	} else {
		sigma2p := indiv.get_varpheno()
		return sigma2p / sigma2env
	}
}

func (cell *Cell) DevCell(G Genome, env Cue) (Cell, error) { //Develops a cell given cue
	var diff float64
	var convindex int

	e_p := NewVec(nenv + ncells) // = env - p0
	g0 := Ones(ngenes)
	f0 := NewVec(ngenes)
	h0 := NewVec(ngenes) //No higher order complexes in embryonic stage
	p0 := NewVec(nenv + ncells)
	ve := NewVec(ngenes)
	vf := NewVec(ngenes)
	vg := NewVec(ngenes)
	vh := NewVec(ngenes)
	vp := NewVec(nenv + ncells)
	f1 := NewVec(ngenes)
	g1 := NewVec(ngenes)
	h1 := NewVec(ngenes)

	convindex = 0
	for nstep := 0; nstep < maxDevStep; nstep++ {
		multMatVec(vf, G.G, g0)
		if withCue { //Model with or without cues
			diffVecs(e_p, env, p0) // env includes noise
			multMatVec(ve, G.E, e_p)
			addVecs(f1, vf, ve)
		} else {
			copy(f1, vf)
		}
		applyFnVec(sigmaf, f1)
		if epig { //Allow or disallow epigenetic layer
			multMatVec(g1, G.F, f1)
			applyFnVec(sigmag, g1)
		} else { //Remove epigenetic layer if false
			copy(g1, f1)
		}
		if hoc { //If layer for higher order complexes is present
			multMatVec(vg, G.H, g1)

			if hoi { //If interactions between higher order complexes is present
				multMatVec(vh, G.J, h0)
				addVecs(h1, vg, vh)
			} else {
				copy(h1, vg)
			}
			applyFnVec(sigmah, h1)
		} else {
			copy(h1, g1) //identity map
		}

		multMatVec_T(vp, G.P, h1)
		applyFnVec(rho, vp)
		diff = DistVecs1(h0, h1) //+ DistVecs1(g0, g1) + DistVecs1(f0, f1) //Stricter convergence criterion requiring convergence in all layers
		//		diff = DistVecs1(p0, vp)
		copy(f0, f1)
		copy(g0, g1)
		copy(h0, h1)
		copy(p0, vp)
		if diff < epsDev { //if criterion is reached
			convindex++ //increment by one
		} else {
			convindex = 0
		}
		if convindex > ccStep {
			cell.NDevStep = nstep //steady state reached
			break
		}
		if nstep == maxDevStep-1 {
			cell.NDevStep = maxDevStep
		}

	}
	//	fmt.Println("PathLen: ", cell.NDevStep)
	//fmt.Println("Phenotype after development:",vpc)
	//fmt.Println("Id after development:",vpid)
	copy(cell.E, env)
	copy(cell.F, f1)
	copy(cell.G, g1)
	copy(cell.H, h1)
	copy(cell.P, vp)

	if cell.NDevStep == maxDevStep {
		return *cell, errors.New("DevCell: did not converge")
	} else {
		return *cell, nil
	}
}

func (body *Body) DevBody(G Genome, envs Cues) (Body, error) {
	var err error = nil
	nenvs := AddNoise2Cues(envs, devNoise)
	for i, cell := range body.Cells {
		cell1, e := cell.DevCell(G, nenvs[i])
		body.Cells[i] = cell1 // Don't forget to update this!
		if e != nil {
			err = e
		}
	}

	return *body, err
}

func (indiv *Indiv) CompareDev(ancenvs, novenvs Cues) Indiv { //Compare developmental process under different conditions

	bodies := indiv.Bodies

	bodies[INoEnv].DevBody(indiv.Genome, ZeroEnvs)
	bodies[IAncEnv].DevBody(indiv.Genome, ancenvs)
	_, errnov := bodies[INovEnv].DevBody(indiv.Genome, novenvs)

	maxdev := 0
	for _, cell := range indiv.Bodies[INovEnv].Cells {
		if cell.NDevStep > maxdev {
			maxdev = cell.NDevStep
		}
	}
	indiv.NDevStep = maxdev

	sse := indiv.get_sse(novenvs)
	if errnov != nil {
		indiv.Fit = 0 //minimum fitness if cells don't converge
	} else {
		indiv.MSE = sse / float64(ncells*(nenv+ncells))
		indiv.Fit = indiv.get_fitness() //Fitness with cues
	}

	// Ignoring convergence/divergence for now
	indiv.AncCuePlas = get_plasticity(indiv.Bodies[INoEnv], indiv.Bodies[IAncEnv])
	indiv.NovCuePlas = get_plasticity(indiv.Bodies[INoEnv], indiv.Bodies[INovEnv])
	indiv.ObsPlas = get_plasticity(indiv.Bodies[IAncEnv], indiv.Bodies[INovEnv])
	indiv.Pp = indiv.get_polyphenism(novenvs)

	return *indiv
}
