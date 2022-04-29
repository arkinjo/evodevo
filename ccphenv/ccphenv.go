package main

// Cross-covariance analysis of state vectors.

import (
	"flag"
	"fmt"
	"log"

	"github.com/arkinjo/evodevo/multicell"
	"gonum.org/v1/gonum/mat"
)

func main() {
	maxpopP := flag.Int("maxpop", 1000, "maximum number of individuals in population")
	ncellsP := flag.Int("ncells", 1, "number of cell types/phenotypes simultaneously trained")

	state0 := flag.String("state0", "P", "State 0 (one of E, F, G, H, P)")
	state1 := flag.String("state1", "P", "State 1 (one of E, F, G, H, P)")
	ienv0 := flag.Int("ienv0", 0, "0 = Ancestral; 1 = Novel environment")
	ienv1 := flag.Int("ienv1", 1, "0 = Ancestral; 1 = Novel environment")

	jsoninP := flag.String("jsonin", "", "json file of input population") //default to empty string

	flag.Parse()

	pop0 := multicell.NewPopulation(*ncellsP, *maxpopP)
	if *jsoninP != "" {
		pop0.FromJSON(*jsoninP)
		multicell.SetParams(pop0.Params)
	} else {
		flag.PrintDefaults()
		log.Fatal("Specify the input JSON file with -jsonin=filename.")
	}

	fmt.Println("# Cross-covariance:", *state0, *ienv0, " vs ", *state1, *ienv1)

	fstates0 := pop0.GetFlatStateVec(*state0, *ienv0)
	fstates1 := pop0.GetFlatStateVec(*state1, *ienv1)
	mstate0, mstate1, ccmat := multicell.GetCrossCov(fstates0, fstates1)

	U, vals, V := runSVD(ccmat)
	fmt.Println("#Singular Values: ", vals)
	project(fstates0, fstates1, mstate0, mstate1, U, V)
	LinearResponse(&pop0)
	PCA_PP(&pop0)
	PCA_PE(&pop0)
}

func PCA_PP(pop *multicell.Population) {
	fp0 := pop.GetFlatStateVec("P", 0)
	fp1 := pop.GetFlatStateVec("P", 1)
	fp := pop.GetFlatStateVec("P", 0)

	for _, p := range fp1 {
		fp = append(fp, p)
	}
	mp, _, cov := multicell.GetCrossCov(fp, fp)
	dim := len(mp)
	mv := mat.NewVecDense(dim, mp)

	evs, U := runEigenSym(cov)

	for i, v := range evs {
		fmt.Printf("PP_val\t%3d\t%e\n", i, v)
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)
	de := mat.NewVecDense(dim, denv)

	for k := range fp0 {
		p0 := mat.NewVecDense(dim, fp0[k])
		p1 := mat.NewVecDense(dim, fp1[k])
		p0.SubVec(p0, mv)
		p1.SubVec(p1, mv)
		fmt.Printf("PP_prj\t%3d", k)
		for i := 0; i < 3; i++ {
			axis := U.ColView(dim - i - 1)
			x := mat.Dot(axis, p0)
			y := mat.Dot(axis, p1)
			if mat.Dot(axis, de) < 0.0 {
				x *= -1.0
				y *= -1.0
			}
			fmt.Printf("\t%e\t%e", x, y)
		}
		fmt.Printf("\n")
	}
}

func PCA_PE(pop *multicell.Population) {
	fp0 := pop.GetFlatStateVec("P", 0)
	fp1 := pop.GetFlatStateVec("P", 1)
	fe0 := pop.GetFlatStateVec("E", 0)
	fe1 := pop.GetFlatStateVec("E", 1)
	fp := pop.GetFlatStateVec("P", 0)
	fe := pop.GetFlatStateVec("E", 0)

	for k, p := range fp1 {
		fp = append(fp, p)
		fe = append(fe, fe1[k])
	}
	pave, eave, cov := multicell.GetCrossCov(fp, fe)
	dim := len(pave)
	mp := mat.NewVecDense(dim, pave)
	me := mat.NewVecDense(dim, eave)

	U, sv, V := runSVD(cov)

	for i, v := range sv {
		fmt.Printf("PE_val\t%3d\t%e\n", i, v)
	}

	env0 := multicell.FlattenEnvs(pop.AncEnvs)
	env1 := multicell.FlattenEnvs(pop.NovEnvs)
	denv := multicell.NewVec(dim)
	multicell.DiffVecs(denv, env1, env0)
	de := mat.NewVecDense(dim, denv)

	for k := range fp0 {
		p0 := mat.NewVecDense(dim, fp0[k])
		p1 := mat.NewVecDense(dim, fp1[k])
		e0 := mat.NewVecDense(dim, fe0[k])
		e1 := mat.NewVecDense(dim, fe1[k])

		p0.SubVec(p0, mp)
		p1.SubVec(p1, mp)

		e0.SubVec(e0, me)
		e1.SubVec(e1, me)

		fmt.Printf("PE_prj\t%3d", k)
		for i := 0; i < 3; i++ {
			u := U.ColView(i)
			v := V.ColView(i)
			xp := mat.Dot(u, p0)
			xe := mat.Dot(v, e0)
			yp := mat.Dot(u, p1)
			ye := mat.Dot(v, e1)
			if mat.Dot(u, de) < 0.0 {
				xp *= -1.0
				yp *= -1.0
			}
			if mat.Dot(v, de) < 0.0 {
				xp *= -1.0
				yp *= -1.0
			}
			fmt.Printf("\t%e\t%e\t%e\t%e", xe, xp, ye, yp)
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
	for i, ci := range pe_cov {
		for j, v := range ci {
			dpp[i] += v * denv[j]
		}
	}

	for i, de := range denv {
		fmt.Printf("LRT\t%2d\t%e\t%e\t%e\n", i, de, dp[i], dpp[i])
	}
}

func project(data0, data1 []multicell.Vec, mean0, mean1 multicell.Vec, U, V *mat.Dense) {
	dim0, _ := U.Dims()
	dim1, _ := V.Dims()

	m0 := mat.NewVecDense(dim0, mean0)
	m1 := mat.NewVecDense(dim1, mean1)

	t0 := mat.NewVecDense(dim0, nil)
	t1 := mat.NewVecDense(dim1, nil)
	for k := range data0 {
		fmt.Printf("Indiv\t%3d", k)
		p0 := mat.NewVecDense(dim0, data0[k])
		for i := 0; i < 3; i++ {
			t0.SubVec(p0, m0)
			axis := U.ColView(i)
			x := mat.Dot(t0, axis)
			fmt.Printf("\t%e", x)
		}
		p1 := mat.NewVecDense(dim1, data1[k])
		for i := 0; i < 3; i++ {
			t1.SubVec(p1, m1)
			axis := V.ColView(i)
			x := mat.Dot(t1, axis)
			fmt.Printf("\t%e", x)
		}
		fmt.Printf("\n")
	}
}

func runSVD(ccmat multicell.Dmat) (*mat.Dense, []float64, *mat.Dense) {
	dim0 := len(ccmat)
	dim1 := len(ccmat[0])
	C := mat.NewDense(dim0, dim1, nil)
	for i, ci := range ccmat {
		for j, v := range ci {
			C.Set(i, j, v)
		}
	}

	var svd mat.SVD
	ok := svd.Factorize(C, mat.SVDFull)
	if !ok {
		log.Fatal("SVD failed.")
	}
	U := mat.NewDense(dim0, dim0, nil)
	V := mat.NewDense(dim1, dim1, nil)
	svd.UTo(U)
	svd.VTo(V)
	vals := svd.Values(nil)

	return U, vals, V
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
