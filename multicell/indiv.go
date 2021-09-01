package multicell

import (
	//"fmt"
	"math"
	"math/rand"
)

type Genome struct { //Genome of an individual
	Ec Spmat //Environment cue effect on epigenome
	Eid Spmat //Environment input of cell type
	F Spmat //Regulatory effect of epigenetic markers on gene expression
	G Spmat //Regulatory effect of gene expression on epigenetic markers
	Hg Spmat //Contribution of gene expression on higher order complexes
	Hh Spmat //Interaction between higher order complexes
	Pc Spmat //Resulting expressed phenotype
	Pid Spmat //Mimiced cell type
	Z Spmat //Gene expression of offspring
}

type Cell struct { //A 'cell' is characterized by its gene expression and phenotype
	E Cue // Environment encountered by cell; id already in cue
	F Vec // Epigenetic markers
	G Vec // Gene expression
	H Vec // Higher order complexes
	P Cue // Phenotype; id already in cue
}

type Cells struct {
	Ctypes	[]Cell // Array of cells of different types
}

type Indiv struct { //An individual as an unicellular organism
	Id       int
	DadId	 int
	MomId	 int
	Genome   Genome
	Copies   []Cells
	Z        Vec // Initial gene expression of offspring
	F0	 	 float64 //Fitness without cues
	Fit		 float64 //Fitness with cues
	Util	 float64 //Fitness Utility of cues 
	CuePlas  float64 //Cue Plasticity
	ObsPlas  float64 //Observed Plasticity
	Pp		 float64 //Degree of polyphenism
}

func NewGenome() Genome { //Generate new genome matrix ensemble
	Ec := NewSpmat(Ngenes, Nenv, GenomeDensity)
	Eid := NewSpmat(Ngenes, ncells, GenomeDensity)
	F := NewSpmat(Ngenes, Ngenes, GenomeDensity)
	G := NewSpmat(Ngenes, Ngenes, GenomeDensity)
	Hg := NewSpmat(Ngenes, Ngenes, GenomeDensity)
	Hh := NewSpmat(Ngenes,Ngenes,GenomeDensity)
	Pc := NewSpmat(Ngenes, Nenv, GenomeDensity)
	Pid := NewSpmat(Ngenes,ncells, GenomeDensity)
	Z := NewSpmat(Ngenes, Ngenes, GenomeDensity)

	genome := Genome{Ec, Eid, F, G, Hg, Hh, Pc, Pid, Z}

	return genome
}

func (parent *Genome) Copy() Genome { //creates copy of parent's genome
	ec := make(Spmat, Ngenes)
	eid := make(Spmat,Ngenes)
	f := make(Spmat, Ngenes)
	g := make(Spmat, Ngenes)
	hg := make(Spmat, Ngenes)
	hh := make(Spmat, Ngenes)
	pc := make(Spmat, Ngenes)
	pid := make(Spmat, Ngenes)
	z := make(Spmat, Ngenes)

	for i, m := range parent.Ec {
		ec[i] = make(map[int]float64)
		for j, v := range m {
			ec[i][j] = v
		}
	}
	for i, m := range parent.Eid {
		eid[i] = make(map[int]float64)
		for j,v := range m {
			eid[i][j] = v
		}
	}
	for i, m := range parent.F {
		f[i] = make(map[int]float64)
		for j, v := range m {
			f[i][j] = v
		}
	}
	for i, m := range parent.G {
		g[i] = make(map[int]float64)
		for j, v := range m {
			g[i][j] = v
		}
	}
	for i, m := range parent.Hg {
		hg[i] = make(map[int]float64)
		for j, v := range m {
			hg[i][j] = v
		}
	}
	for i, m := range parent.Hh {
		hh[i] = make(map[int]float64)
		for j, v := range m {
			hh[i][j] = v
		}
	}
	for i, m := range parent.Pc {
		pc[i] = make(map[int]float64)
		for j, v := range m {
			pc[i][j] = v
		}
	}
	for i, m := range parent.Pid {
		pid[i] = make(map[int]float64)
		for j, v := range m {
			pid[i][j] = v
		}
	}
	for i, m := range parent.Z {
		z[i] = make(map[int]float64)
		for j, v := range m {
			z[i][j] = v
		}
	}


	genome := Genome{ec,eid,f,g,hg,hh,pc,pid,z}

	return genome
}

func DiffGenomes(Gout, G1, G0 Genome){ //Elementwise difference between two genomes
	for i := 0 ; i < Ngenes ; i++ {
		for j := 0; j < Nenv; j++ {
			Gout.Ec[i][j] = G1.Ec[i][j] - G0.Ec[i][j]
		}
		for j := 0; j < ncells; j++ {
			Gout.Eid[i][j] = G1.Eid[i][j] - G0.Eid[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.F[i][j] = G1.F[i][j] - G0.F[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.G[i][j] = G1.G[i][j] - G0.G[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.Hg[i][j] = G1.Hg[i][j] - G0.Hg[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.Hh[i][j] = G1.Hh[i][j] - G0.Hh[i][j]
		}
		for j := 0; j < Nenv; j++ {
			Gout.Pc[i][j] = G1.Pc[i][j] - G0.Pc[i][j]
		}
		for j := 0; j < ncells; j++ {
			Gout.Pid[i][j] = G1.Pid[i][j] - G1.Pid[i][j]
		}
		for j := 0; j < Ngenes; j++ {
			Gout.Z[i][j] = G1.Z[i][j] - G0.Z[i][j]
		}
	}
}

func(G *Genome) NormalizeGenome() Genome{
	lambda2 := 0.0
	eG := G.Copy()
	for _,m := range G.Ec{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Eid{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.F{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.G{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Hg{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Hh{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Pc{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Pid{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	for _,m := range G.Z{
		for _,v := range m{
			lambda2 += v*v
		}
	}
	lambda := math.Sqrt(lambda2)
	for i,m := range eG.Ec{
		for j := range m{
			eG.Ec[i][j] = eG.Ec[i][j]/lambda
		}
	}
	for i,m := range eG.Eid{
		for j := range m{
			eG.Eid[i][j] = eG.Eid[i][j]/lambda
		}
	}
	for i,m := range eG.F{
		for j := range m{
			eG.F[i][j] = eG.F[i][j]/lambda
		}
	}
	for i,m := range eG.G{
		for j := range m{
			eG.G[i][j] = eG.G[i][j]/lambda
		}
	}
	for i,m := range eG.Hg{
		for j := range m{
			eG.Hg[i][j] = eG.Hg[i][j]/lambda
		}
	}
	for i,m := range eG.Hh{
		for j := range m{
			eG.Hh[i][j] = eG.Hh[i][j]/lambda
		}
	}
	for i,m := range eG.Pc{
		for j := range m{
			eG.Pc[i][j] = eG.Pc[i][j]/lambda
		}
	}
	for i,m := range eG.Pid{
		for j := range m{
			eG.Pid[i][j] = eG.Pid[i][j]/lambda
		}
	}
	for i,m := range eG.Z{
		for j := range m{
			eG.Z[i][j] = eG.Z[i][j]/lambda
		}
	}
	return eG
}

func NewCell(id int) Cell { //Creates a new cell
	e := NewCue(Nenv,id)
	f := NewVec(Ngenes)
	g := NewVec(Ngenes)
	h := NewVec(Ngenes)
	p := NewCue(Nenv,id)

	cell := Cell{e,f,g,h,p}

	return cell
}

func NewCells(ncells int) Cells { // Creates an array of new cells of length Ncells
	Clist := make([]Cell, ncells)
	for id := range Clist {
		Clist[id] = NewCell(id) //Initialize each cell
	}
	Cellarray := Cells{Clist}
	return Cellarray
}

func NewIndiv(id int) Indiv { //Creates a new individual
	genome := NewGenome()
	cellcopies := make([]Cells,3)
	for i := range(cellcopies) {
		cellcopies[i] = NewCells(ncells)
	}

	z := NewVec(Ngenes)

	indiv := Indiv{id, 0, 0, genome, cellcopies, z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	return indiv
}

func (indiv *Indiv) Mutate() { //Mutates genome of an individual
	ipos := rand.Intn(GeneLength)
	if ipos < Nenv {
		mutateSpmat(indiv.Genome.Ec,Nenv)
	} else if ipos < Nenv + ncells {
		mutateSpmat(indiv.Genome.Eid,ncells)
	} else if ipos < Nenv + Ngenes + ncells {
		mutateSpmat(indiv.Genome.F,Ngenes)
	} else if ipos < Nenv + 2*Ngenes + ncells{
		mutateSpmat(indiv.Genome.G,Ngenes)
	} else if ipos < Nenv + 3*Ngenes + ncells{
		mutateSpmat(indiv.Genome.Hg,Ngenes)
	} else if ipos < Nenv + 4*Ngenes + ncells{
		mutateSpmat(indiv.Genome.Hh,Ngenes)
	} else if ipos < 2*Nenv + 4*Ngenes + ncells{
		mutateSpmat(indiv.Genome.Pc,Nenv)
	} else if ipos < 2*Nenv + 4*Ngenes + 2*ncells {
		mutateSpmat(indiv.Genome.Eid,ncells)
	} else {
		mutateSpmat(indiv.Genome.Z,Ngenes)
	}

	//May need variations based on model used (E.g. no mutations to E when withCue == false)

	return
}

func Crossover(dadg, momg *Genome, dadz , momz Vec) (Genome, Genome, Vec, Vec) { //Crossover
	ng0 := dadg.Copy()
	ng1 := momg.Copy()
	nz0 := dadz
	nz1 := momz

	for i := 0; i < Ngenes; i++ {
		r := rand.Float64()
		if r < 0.5 {
			ec := ng0.Ec[i]
			eid := ng0.Eid[i]
			f := ng0.F[i]
			g := ng0.G[i]
			hg := ng0.Hg[i]
			hh := ng0.Hh[i]
			pc := ng0.Pc[i]
			pid := ng0.Pid[i]
			z := ng0.Z[i]
			z0 := nz0[i]

			ng0.Ec[i] = ng1.Ec[i]
			ng0.Eid[i] = ng1.Eid[i]
			ng0.F[i] = ng1.F[i]
			ng0.G[i] = ng1.G[i]
			ng0.Hg[i] = ng1.Hg[i]
			ng0.Hh[i] = ng1.Hh[i]
			ng0.Pc[i] = ng1.Pc[i]
			ng0.Pid[i] = ng1.Pid[i]
			ng0.Z[i] = ng1.Z[i]
			nz0[i] = nz1[i]

			ng1.Ec[i] = ec
			ng1.Eid[i] = eid
			ng1.F[i] = f
			ng1.G[i] = g
			ng1.Hg[i] = hg
			ng1.Hh[i] = hh
			ng1.Pc[i] = pc
			ng1.Pid[i] = pid
			ng1.Z[i] = z
			nz1[i] = z0
		}
	}

	return ng0, ng1, nz0, nz1
}

func Mate(dad, mom *Indiv) (Indiv, Indiv) { //Generates offspring
	genome0, genome1, g0, g1 :=
		Crossover(&dad.Genome, &mom.Genome, dad.Z, mom.Z)

	cells0 := make([]Cells,3)
	for i := range cells0 {
		cells0[i] = NewCells(ncells)
	}
	cells1 := make([]Cells,3)
	for i := range cells1 {
		cells1[i] = NewCells(ncells)
	}

	kid0 := Indiv{dad.Id, dad.Id, mom.Id, genome0,  cells0, g0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	kid1 := Indiv{mom.Id, dad.Id, mom.Id, genome1,  cells1, g1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	kid0.Mutate()
	kid1.Mutate()

	return kid0, kid1
}

/*
func (cell *Cell) get_fitness(env Cue) float64 {
	d2 := dist2Vecs(env.C,cell.P.C)
	return math.Exp(-s * d2/float64(Ncells))
}
*/

func (cells *Cells) get_fitness(envs Cues) float64 {
	d2 := 0.0
	env := make([]float64,Nenv)
	p := make([]float64,Nenv)
	id := make([]float64,ncells)
	idp := make([]float64,ncells)
	for i := range cells.Ctypes{
		
		env = envs.Es[i].C
		p = cells.Ctypes[i].P.CopyCue().C
		d2 += dist2Vecs(p,env)

		//Also add difference of id of cell type
		id = envs.Es[i].id
		idp = cells.Ctypes[i].P.CopyCue().id
		//fmt.Println("env id:",id)
		//fmt.Println("cell id:",idp)
		d2 += dist2Vecs(idp,id)
	}
	// s := ScaleSelStrength(Ncell)
	// return(...(-s*d2))
	return math.Exp(-selStrength * d2)
}


func (indiv *Indiv) get_cue_plasticity() float64 { //cue plasticity of individual
	d2 := 0.0
	p0 := NewVec(Nenv)
	p := NewVec(Nenv)
	Copies := indiv.Copies
	for i,cell := range Copies[2].Ctypes {
		p = cell.P.C
		p0 = Copies[0].Ctypes[i].P.C
		d2 += dist2Vecs(p,p0)
	}
	d2 = dist2Vecs(p,p0)/float64(Nenv*ncells) //Divide by number of phenotypes to normalize
	return d2
}

func (indiv *Indiv) get_obs_plasticity() float64 { //cue plasticity of individual
	d2 := 0.0
	p0 := NewVec(Nenv)
	p := NewVec(Nenv)
	Copies := indiv.Copies
	for i,cell := range Copies[2].Ctypes {
		p = cell.P.C
		p0 = Copies[1].Ctypes[i].P.C
		d2 += dist2Vecs(p,p0)
	}
	d2 = dist2Vecs(p,p0)/float64(Nenv*ncells) //Divide by number of phenotypes to normalize
	return d2
}

func(indiv *Indiv) get_vp() float64 { //Get sum of elementwise variance of phenotype
	pvec := make([]Cue,0)
	for _,c := range indiv.Copies[2].Ctypes {
		pvec = append(pvec,c.P) //Note: To be used AFTER development
	}
	phens := Cues{pvec}
	sigma2p := phens.GetCueVar()
	return sigma2p
}

func(indiv *Indiv) get_pp(envs Cues) float64 { //Degree of polyphenism of individual; normalize with variance of environment cue
	sigma2env := envs.GetCueVar()
	if sigma2env == 0 {
		return 0
	} else {
		sigma2p := indiv.get_vp()
		return sigma2p/sigma2env 
	}
}

func (cell *Cell) DevCell(G Genome, g0 Vec, env Cue) Cell { //Develops a cell given cue
	var diff float64
	
	h0 := make([]float64,Ngenes) //No higher order complexes in embryonic stage
	vec := make([]float64,Ngenes)
	veid := make([]float64,Ngenes)
	vf := make([]float64,Ngenes)
	vg := g0
	vh := make([]float64,Ngenes)
	vpc := make([]float64,Nenv)
	vpid := make([]float64,ncells)
	e1 := make([]float64,Ngenes)
	f1 := make([]float64,Ngenes)
	g1 := make([]float64,Ngenes)
	h1 := make([]float64,Ngenes)

	for nstep := 0 ; nstep < MaxDevStep ; nstep++ {
		multMatVec(vec,G.Ec,cell.E.C)
		multMatVec(veid,G.Eid,cell.E.id)
		addVecs(e1,vec,veid)
		multMatVec(vf,G.G,g0)
		if WithCue { //Model with or without cues
			addVecs(f1,vf,e1)
		} else {
			f1 = vf
		}
		applyFnVec(sigmaf,f1)
		if Epig { //Allow or disallow epigenetic layer
			multMatVec(g1,G.F,f1)
			applyFnVec(sigmag,g1)
		} else { //Remove epigenetic layer if false
			g1 = f1
		}
		if HOC { //If layer for higher order complexes is present
			multMatVec(vg,G.Hg,g1)
			multMatVec(vh,G.Hh,h0)
			if HOI { //If interactions between higher order complexes is present
				addVecs(h1,vg,vh)
			} else {
				h1 = vg
			}
		} else {
			h1 = g1
		}
		applyFnVec(sigmah,h1)
		multMatVec_T(vpc,G.Pc,h1)
		applyFnVec(rho,vpc)
		multMatVec_T(vpid,G.Pid,h1)
		applyFnVec(rho,vpid)
		diff = dist2Vecs(vg,g0)
		g0 = g1
		h0 = h1
		if diff < epsDev{
			break
		}

	}
	//fmt.Println("Phenotype after development:",vpc)
	//fmt.Println("Id after development:",vpid)


	cell.E = env
	cell.F = f1
	cell.G = g1
	cell.H = h1
	cell.P = Cue{vpid,vpc}
	
	return *cell
}

func (cells *Cells) DevCells(G Genome, g0 Vec, envs Cues) Cells {
	env := NewCue(Nenv,0)
	
	for i := range cells.Ctypes {
		env = envs.Es[i]
		cells.Ctypes[i].DevCell(G, g0, env)
	}

	return *cells
}

func (indiv *Indiv) CompareDev(env, env0 Cues) Indiv { //Compare developmental process under different conditions

	devenv := env.AddNoise(DevNoise)
	devenv0 := env0.AddNoise(DevNoise)
	selenv := env.AddNoise(EnvNoise)
	Clist := indiv.Copies
	ncells := len(Clist[2].Ctypes)
	nenv := len(Clist[2].Ctypes[0].E.C)
	zero := NewCues(ncells,nenv)

	Clist[0].DevCells(indiv.Genome, indiv.Z, zero) //Develop without cues
	Clist[1].DevCells(indiv.Genome, indiv.Z, devenv0) //Develop in ancestral (previous) environment
	Clist[2].DevCells(indiv.Genome, indiv.Z, devenv) //Develop in novel (present) environment
	
	multMatVec_T(indiv.Z,indiv.Genome.Z,Clist[2].Ctypes[0].G)
	indiv.F0 = Clist[0].get_fitness(selenv)  //Fitness without cues
	indiv.Fit = Clist[2].get_fitness(selenv) //Fitness with cues
	indiv.Util = indiv.Fit-indiv.F0

	indiv.CuePlas = indiv.get_cue_plasticity()
	indiv.ObsPlas = indiv.get_obs_plasticity()
	indiv.Pp = indiv.get_pp(devenv)

	return *indiv
}
