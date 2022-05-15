package multicell

import (
	//	"errors"
	//	"fmt"
	"log"
	"math"
	//	"math/rand"
)

type Cell struct { //A 'cell' is characterized by its gene expression and phenotype
	E        Vec     // Env + noise
	F        Vec     // Epigenetic markers
	G        Vec     // Gene expression
	H        Vec     // Higher order complexes
	P, Pvar  Vec     // moving average and variance of P (P is already EMA)
	PErr     float64 // ||e - p||_1
	NDevStep int     // Developmental path length

}

type Body struct { //Do we want to reimplement this?
	Genome   Genome
	Cells    []Cell // Array of cells of different types
	PErr     float64
	NDevStep int
}

const (
	IAncEnv = iota // Previous env
	INovEnv        // Current env
)

const ( // Index for cell state vectors
	CellE = iota
	CellF
	CellG
	CellH
	CellP
)

type Indiv struct { //An individual as an unicellular organism
	Id         int
	DadId      int
	MomId      int
	Bodies     []Body  //IAncEnv, INovEnv (see above const.)
	Fit        float64 //Fitness with cues
	WagFit     float64 //Wagner relative fitness
	Plasticity float64 //Observed Plasticity
	Dp1e1      float64 // ||p(e1) - e1||
	Dp0e0      float64 // ||p(e0) - e0||
	Dp1e0      float64 // ||p(e1) - e0||
	Dp0e1      float64 // ||p(e0) - e1||
}

func NewCell(id int) Cell { //Creates a new cell given id of cell.
	e := NewVec(nenv + ncells)
	f := NewVec(ngenes)
	g := NewVec(ngenes)
	h := NewVec(ngenes)
	p := NewCue(nenv, id)
	pv := NewVec(nenv + ncells)
	cell := Cell{e, f, g, h, p, pv, 0.0, 0}

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
	cell1.Pvar = CopyVec(cell.Pvar)
	cell1.PErr = cell.PErr
	cell1.NDevStep = cell.NDevStep

	return cell1
}

func (cell *Cell) GetState(ivec string) Vec {
	switch ivec {
	case "E":
		return cell.E
	case "F":
		return cell.F
	case "G":
		return cell.G
	case "H":
		return cell.H
	case "P":
		return cell.P
	default:
		log.Fatal("Cell.GetState: Unknown state vector")
	}

	return nil // neven happens
}

func NewBody(ncells int) Body {
	genome := NewGenome()
	cells := make([]Cell, ncells)
	for id := range cells {
		cells[id] = NewCell(id) //Initialize each cell
	}
	return Body{genome, cells, 0, 0}
}

func (body *Body) Copy() Body {
	body1 := NewBody(ncells)
	body1.PErr = body.PErr
	body1.NDevStep = body.NDevStep
	body1.Genome = body.Genome.Copy()
	for i, cell := range body.Cells {
		body1.Cells[i] = cell.Copy()
	}
	return body1
}

func NewIndiv(id int) Indiv { //Creates a new individual
	bodies := make([]Body, 2)
	for i := range bodies {
		bodies[i] = NewBody(ncells)
	}

	indiv := Indiv{id, 0, 0, bodies, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	return indiv
}

func (indiv *Indiv) Copy() Indiv { //Deep copier
	indiv1 := NewIndiv(indiv.Id)
	indiv1.DadId = indiv.DadId
	indiv1.MomId = indiv.MomId
	for i, body := range indiv.Bodies {
		indiv1.Bodies[i] = body.Copy()
	}
	indiv1.Fit = indiv.Fit
	indiv1.Plasticity = indiv.Plasticity
	indiv1.Dp1e1 = indiv.Dp1e1
	indiv1.Dp0e0 = indiv.Dp0e0
	indiv1.Dp1e0 = indiv.Dp1e0
	indiv1.Dp0e1 = indiv.Dp0e1

	return indiv1
}

func Mate(dad, mom *Indiv) (Indiv, Indiv) { //Generates offspring
	genome0 := dad.Bodies[INovEnv].Genome.Copy()
	genome1 := mom.Bodies[INovEnv].Genome.Copy()
	CrossoverSpmats(genome0.E, genome1.E)
	CrossoverSpmats(genome0.F, genome1.F)
	CrossoverSpmats(genome0.G, genome1.G)
	CrossoverSpmats(genome0.H, genome1.H)
	CrossoverSpmats(genome0.J, genome1.J)
	CrossoverSpmats(genome0.P, genome1.P)

	bodies0 := make([]Body, 3)
	for i := range bodies0 {
		bodies0[i] = NewBody(ncells)
		bodies0[i].Genome = genome0
	}
	bodies1 := make([]Body, 3)
	for i := range bodies1 {
		bodies1[i] = NewBody(ncells)
		bodies1[i].Genome = genome1
	}

	kid0 := Indiv{dad.Id, dad.Id, mom.Id, bodies0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	kid1 := Indiv{mom.Id, dad.Id, mom.Id, bodies1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	kid0.Bodies[INovEnv].Genome.Mutate()
	kid1.Bodies[INovEnv].Genome.Mutate()

	return kid0, kid1
}

func (indiv *Indiv) getPErr(ienv int) float64 {
	return indiv.Bodies[ienv].PErr
}

func (indiv *Indiv) getNDevStep(ienv int) int {
	return indiv.Bodies[ienv].NDevStep
}

func (indiv *Indiv) getFitness() float64 { //fitness in novel/present environment
	ndevstep := indiv.getNDevStep(INovEnv)
	if ndevstep == maxDevStep {
		return 0.0
	}

	fdev := float64(ndevstep) / selDevStep
	ferr := indiv.getPErr(INovEnv) * baseSelStrength
	rawfit := math.Exp(-(ferr + fdev))
	return rawfit
}

func getPlasticity(body0, body1 Body) float64 { //cue plasticity of individual
	d2 := 0.0
	for i, cell := range body0.Cells {
		d2 += Dist2Vecs(cell.P, body1.Cells[i].P)
	}

	return d2 / float64(ncells*(ncells+nenv))
}

func getPEDiff(body Body, envs Cues) float64 {
	diff := 0.0
	for i, c := range body.Cells {
		diff += DistVecs1(c.P, envs[i])
	}
	return diff / float64(ncells*(ncells+nenv))
}

func (indiv *Indiv) getVarpheno() float64 { //Get sum of elementwise variance of phenotype
	pvec := make([]Cue, 0)
	for _, c := range indiv.Bodies[INovEnv].Cells {
		pvec = append(pvec, c.P) //Note: To be used AFTER development
	}
	sigma2p := GetCueVar(pvec)
	return sigma2p
}

func (cell *Cell) getPscale() Cue {
	vt := NewVec(nenv + ncells)
	evar := devNoise * devNoise
	for i, v := range cell.Pvar {
		vt[i] = v / (v + evar)
	}
	log.Println("Kalman gain", vt)
	return vt
}

func (cell *Cell) updatePEMA(pnew Vec) {
	for i, pi := range pnew {
		d := pi - cell.P[i]
		incr := alphaEMA * d
		cell.P[i] += incr
		cell.Pvar[i] = (1 - alphaEMA) * (cell.Pvar[i] + d*incr)
	}
}

func (cell *Cell) DevCell(G Genome, env Cue) Cell { //Develops a cell given cue
	e_p := NewVec(nenv + ncells) // = env - p0
	g0 := Ones(ngenes)
	f0 := NewVec(ngenes)
	h0 := NewVec(ngenes) //No higher order complexes in embryonic stage
	ve := NewVec(ngenes)
	vf := NewVec(ngenes)
	vg := NewVec(ngenes)
	vh := NewVec(ngenes)
	p1 := NewVec(nenv + ncells)
	f1 := NewVec(ngenes)
	g1 := NewVec(ngenes)
	h1 := NewVec(ngenes)

	AddNoise2Cue(cell.E, env, devNoise)

	for nstep := 1; nstep <= maxDevStep; nstep++ {
		MultMatVec(vf, G.G, g0)
		if withE { //Model with or without cues
			if pheno_feedback { //If feedback is allowed

				DiffVecs(e_p, cell.E, cell.P)

				// Kalman gain
				//pscale := cell.getPscale()
				//multVecVec(e_p, pscale, e_p)

				MultMatVec(ve, G.E, e_p)
			} else {
				MultMatVec(ve, G.E, cell.E)
			}
			AddVecs(f1, vf, ve)
		} else {
			copy(f1, vf)
		}
		applyFnVec(sigmaf, f1)
		if withF { //Allow or disallow epigenetic layer
			MultMatVec(g1, G.F, f1)
			applyFnVec(sigmag, g1)
		} else { //Remove epigenetic layer if false
			copy(g1, f1)
		}
		if withH { //If layer for higher order complexes is present
			MultMatVec(vg, G.H, g1)

			if withJ { //If interactions between higher order complexes is present
				MultMatVec(vh, G.J, h0)
				AddVecs(h1, vg, vh)
			} else {
				copy(h1, vg)
			}
			applyFnVec(sigmah, h1)
		} else {
			copy(h1, g1) //identity map
		}

		MultMatVec(p1, G.P, h1)
		applyFnVec(rho, p1)
		copy(f0, f1)
		copy(g0, g1)
		copy(h0, h1)
		cell.updatePEMA(p1)
		diff := 0.0
		for _, v := range cell.Pvar {
			diff += v
		}
		diff /= float64(len(cell.Pvar))
		cell.NDevStep = nstep
		if diff < epsDev {
			break
		}
	}
	copy(cell.F, f1)
	copy(cell.G, g1)
	copy(cell.H, h1)
	cell.PErr = DistVecs1(cell.P, env) / cueMag
	//cell.PErr = Dist2Vecs(cell.P, env) / (cueMag * cueMag)

	return *cell
}

func (body *Body) DevBody(envs Cues) Body {
	sse := 0.0
	maxdev := 0

	for i, cell := range body.Cells {
		body.Cells[i] = cell.DevCell(body.Genome, envs[i])
		sse += cell.PErr
		if cell.NDevStep > maxdev {
			maxdev = cell.NDevStep
		}
		if cell.NDevStep > maxDevStep {
			log.Println("NDevStep greater than limit: ", cell.NDevStep)
		}
	}

	body.PErr = sse / float64(ncells*(ncells+nenv))
	body.NDevStep = maxdev

	return *body
}

func (indiv *Indiv) Develop(ancenvs, novenvs Cues) Indiv { //Compare developmental process under different conditions
	indiv.Bodies[IAncEnv].DevBody(ancenvs)
	indiv.Bodies[INovEnv].DevBody(novenvs)

	indiv.Fit = indiv.getFitness()

	indiv.Plasticity = getPlasticity(indiv.Bodies[IAncEnv], indiv.Bodies[INovEnv])
	indiv.Dp1e1 = getPEDiff(indiv.Bodies[INovEnv], novenvs)
	indiv.Dp0e0 = getPEDiff(indiv.Bodies[IAncEnv], ancenvs)
	indiv.Dp1e0 = getPEDiff(indiv.Bodies[INovEnv], ancenvs)
	indiv.Dp0e1 = getPEDiff(indiv.Bodies[IAncEnv], novenvs)
	return *indiv
}
