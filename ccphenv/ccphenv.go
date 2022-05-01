package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")

	jsoninP := flag.String("jsonin", "", "json file of input population") //default to empty string

	flag.Parse()

	pop := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsoninP != "" {
		pop.FromJSON(*jsoninP)
		multicell.SetParams(pop.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	dim := len(env0)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)

	e0 := pop.GetFlatStateVec("E", 0)
	e1 := pop.GetFlatStateVec("E", 1)

	p0 := pop.GetFlatStateVec("P", 0)
	p1 := pop.GetFlatStateVec("P", 1)
	mp0 := multicell.GetMeanVec(p0)
	mp1 := multicell.GetMeanVec(p1)
	dp := multicell.NewVec(dim)
	multicell.DiffVecs(dp, mp1, mp0)

	deltaP := make([]multicell.Vec, 0)
	deltaE := make([]multicell.Vec, 0)
	for k, p := range p1 {
		d := multicell.NewVec(dim)
		multicell.DiffVecs(d, p, p0[k])
		deltaP = append(deltaP, d)

		e := multicell.NewVec(dim)
		multicell.DiffVecs(e, e1[k], e0[k])
		deltaE = append(deltaE, e)
	}

	dotPE := multicell.DotVecs(denv, dp)
	dotPP := multicell.DotVecs(dp, dp)
	pefn2 := 0.0
	ppfn2 := 0.0
	for _, p := range dp {
		for _, p2 := range dp {
			ppfn2 += (p * p2) * (p * p2)
		}
		for _, e := range denv {
			pefn2 += (p * e) * (p * e)
		}
	}
	dotEE := multicell.DotVecs(denv, denv)
	eefn2 := 0.0
	for _, e := range denv {
		for _, e2 := range denv {
			eefn2 += (e * e2) * (e * e2)
		}
	}

	dirE := mat.NewVecDense(dim, denv)
	dirE.ScaleVec(1.0/dirE.Norm(2), dirE)
	dirP := mat.NewVecDense(dim, dp)
	dirP.ScaleVec(1.0/dirP.Norm(2), dirP)

	deltaXX("<dPdP>", dirE, deltaP, deltaP)
	fmt.Printf("<dP><dP>_FN2,Tr\t%e\t%e\n", ppfn2, dotPP)
	Project("<dDPdDP>", dirE, dirP, deltaP, deltaP)
	Project("<DP0DP0>", dirE, dirP, p0, p0)
	Project("<DP0DP1>", dirE, dirP, p0, p1)
	Project("<DP1DP0>", dirE, dirP, p1, p0)
	Project("<DP1DP1>", dirE, dirP, p1, p1)

	deltaXX("<dPdE>", dirE, deltaP, deltaE)
	fmt.Printf("<dP><dE>FN2,Tr\t%e\t%e\n", pefn2, dotPE)
	Project("<dDPdDE>", dirE, dirP, deltaP, deltaE)
	Project("<DP0DE0>", dirE, dirP, p0, e0)
	Project("<DP0DE1>", dirE, dirP, p0, e1)
	Project("<DP1DE0>", dirE, dirP, p1, e0)
	Project("<DP1DE1>", dirE, dirP, p1, e1)

	deltaXX("<dEdE>", dirE, deltaE, deltaE)
	fmt.Printf("<dE><dE>FN2,Tr\t%e\t%e\n", eefn2, dotEE)
	Project("<dDEdDE>", dirE, dirP, deltaE, deltaE)
	Project("<DE0DE0>", dirE, dirP, e0, e0)
	Project("<DE0DE1>", dirE, dirP, e0, e1)
	Project("<DE1DE0>", dirE, dirP, e1, e0)
	Project("<DE1DE1>", dirE, dirP, e1, e1)

	mixp := make([]multicell.Vec, 0)
	mixe := make([]multicell.Vec, 0)
	for k := range p0 {
		mixp = append(mixp, p0[k])
		mixe = append(mixe, e0[k])
	}
	for k := range p1 {
		mixp = append(mixp, p1[k])
		mixe = append(mixe, e1[k])
	}
	Project("mixPP", dirE, dirP, mixp, mixp)
	Project("mixPE", dirE, dirP, mixp, mixe)

	LinearResponse(&pop)
}

func deltaXX(label string, dir *mat.VecDense, data0, data1 [][]float64) {
	dim0 := len(data0[0])
	dim1 := len(data1[0])
	cov := multicell.NewDmat(dim0, dim1)
	for k, p := range data0 {
		for i, v0 := range p {
			for j, v1 := range data1[k] {
				cov[i][j] += v0 * v1
			}
		}
	}
	fn := 1.0 / float64(len(data0))
	for i, ci := range cov {
		for j := range ci {
			cov[i][j] *= fn
		}
	}

	trace := 0.0
	for i, v := range cov {
		trace += v[i]
	}

	U, vals, V := multicell.GetSVD(cov)
	fn2 := 0.0
	for _, v := range vals {
		fn2 += v * v
	}
	fmt.Printf("%s_FN2,Tr\t%e\t%e\n", label, fn2, trace)

	for i, v := range vals {
		fmt.Printf("%s_val\t%3d\t%e\n", label, i, v)
	}

	for k := range data0 {
		p0 := mat.NewVecDense(dim0, data0[k])
		p1 := mat.NewVecDense(dim1, data1[k])
		fmt.Printf("%s_prj\t%3d", label, k)
		for i := 0; i < 3; i++ {
			u := U.ColView(i)
			v := V.ColView(i)
			y := mat.Dot(u, p0)
			x := mat.Dot(v, p1)
			if mat.Dot(u, dir) < 0.0 {
				x *= -1.0
				y *= -1.0
			}
			fmt.Printf("\t%e\t%e", x, y)
		}
		fmt.Printf("\n")
	}
}

func LinearResponse(pop *multicell.Population) {
	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	dim := len(env0)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)

	fe0 := pop.GetFlatStateVec("E", 0)
	fp0 := pop.GetFlatStateVec("P", 0)
	fp1 := pop.GetFlatStateVec("P", 1)
	p0_ave, _, pe_cov := multicell.GetCrossCov(fp0, fe0)
	p1_ave := multicell.GetMeanVec(fp1)
	dp := multicell.NewVec(dim)
	multicell.DiffVecs(dp, p1_ave, p0_ave)

	dpp := multicell.NewVec(dim)
	sigma := pop.Params.SDNoise
	sigma2 := (sigma * sigma)
	for i, ci := range pe_cov {
		for j, v := range ci {
			dpp[i] += v * denv[j]
		}
	}

	for i, de := range denv {
		fmt.Printf("LRT\t%2d\t%e\t%e\t%e\n", i, de, dp[i], dpp[i]/sigma2)
	}
}

func Project(label string, dirE, dirP *mat.VecDense, data0, data1 [][]float64) {

	mean0, mean1, ccmat := multicell.GetCrossCov(data0, data1)

	trace := 0.0
	for i, v := range ccmat {
		trace += v[i]
	}

	U, vals, V := multicell.GetSVD(ccmat)

	fnorm2 := 0.0
	for _, v := range vals {
		fnorm2 += v * v
	}
	fmt.Printf("%s_FN2,Tr\t%e\t%e\n", label, fnorm2, trace)
	dim0, _ := U.Dims()
	dim1, _ := V.Dims()

	m0 := mat.NewVecDense(dim0, mean0)
	m1 := mat.NewVecDense(dim1, mean1)

	t0 := mat.NewVecDense(dim0, nil)
	t1 := mat.NewVecDense(dim1, nil)
	for i, v := range vals {
		fmt.Printf("%s_vals\t%d\t%e\n", label, i, v)
	}

	fmt.Println("#           \tcomp\tu.e     \tu.p     \tv.e     \tv.p")
	for i := 0; i < dim0; i++ {
		u := U.ColView(i)
		v := V.ColView(i)
		ue := math.Abs(mat.Dot(dirE, u))
		up := math.Abs(mat.Dot(dirP, u))
		ve := math.Abs(mat.Dot(dirE, v))
		vp := math.Abs(mat.Dot(dirP, v))
		fmt.Printf("%s_aliUV_EP\t%d\t%e\t%e\t%e\t%e\n", label, i, ue, up, ve, vp)
	}

	fmt.Printf("#<YX>   \tcomp")
	for i := 0; i < 3; i++ {
		fmt.Printf("\tX.u%-5d\tY.v%-5d", i, i)
	}
	fmt.Printf("\n")
	for k := range data0 {
		fmt.Printf("%s_prj\t%3d", label, k)
		p0 := mat.NewVecDense(dim0, data0[k])
		p1 := mat.NewVecDense(dim1, data1[k])
		for i := 0; i < 3; i++ {
			t0.SubVec(p0, m0)
			u := U.ColView(i)
			y := mat.Dot(t0, u)
			t1.SubVec(p1, m1)
			v := V.ColView(i)
			x := mat.Dot(t1, v)
			if mat.Dot(dirE, u) < 0.0 {
				x *= -1
				y *= -1
			}
			fmt.Printf("\t%e \t%e ", x, y)
		}
		fmt.Printf("\n")
	}
}

func runEigenSym(cov multicell.Dmat) ([]float64, *mat.Dense) {
	dim := len(cov)
	a := mat.NewSymDense(dim, nil)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			a.SetSym(i, j, cov[i][j])
		}
	}

	var eig mat.EigenSym
	ok := eig.Factorize(a, true)
	if !ok {
		log.Fatal("EigenSym failed.")
	}

	U := mat.NewDense(dim, dim, nil)
	eig.VectorsTo(U)

	return eig.Values(nil), U
}

func getTest3by3Matrix() multicell.Dmat {

	dmat := multicell.NewDmat(3, 3)
	dmat[0][0] = 1
	dmat[0][1] = 2
	dmat[0][2] = -3
	dmat[1][0] = 2
	dmat[1][1] = 3
	dmat[1][2] = 4
	dmat[2][0] = -3
	dmat[2][1] = 4
	dmat[2][2] = 5

	return dmat
}
