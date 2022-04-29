package multicell

import (
	"fmt"
	"log"
	"os"
)

//Dumps genealogy of population for an epoch into a dot file, going backwards in time. Returns number of reproducing population
func DOT_Genealogy(dotfilename, popfilename string, ngen, npop int) []int {
	var npars int
	var id, dadid, momid string

	nanctraj := []int{}
	rnanctraj := []int{}
	pop := NewPopulation(ncells, npop)
	genfile := fmt.Sprintf("../analysis/%s.dot", dotfilename)
	fout, err := os.OpenFile(genfile, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout, "digraph G {")
	kids := make(map[int]bool)
	for i := 0; i < npop; i++ {
		kids[i] = true
	}
	pars := make(map[int]bool)
	for gen := ngen; gen > 0; gen-- {
		jfilename := fmt.Sprintf("../pops/%s_%3.3d.json", popfilename, gen)
		pop.FromJSON(jfilename)

		for _, indiv := range pop.Indivs {
			pars[indiv.DadId] = true
			pars[indiv.MomId] = true
			id = fmt.Sprintf("g%d:id%d", pop.Gen, indiv.Id)
			dadid = fmt.Sprintf("g%d:id%d", pop.Gen-1, indiv.DadId)
			momid = fmt.Sprintf("g%d:id%d", pop.Gen-1, indiv.MomId)
			fmt.Fprintf(fout, "\t \"%s\"-> {\"%s\", \"%s\"}\n", id, dadid, momid)
		}
		npars = len(pars)
		rnanctraj = append(rnanctraj, npars)
		kids = make(map[int]bool)
		for i := range pars {
			kids[i] = pars[i]
		}
		pars = make(map[int]bool) //re-initialize
	}
	fmt.Fprintln(fout, "}")
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}

	//Order of proptraj is recorded backwards, need to reverse it
	copy(nanctraj, rnanctraj)
	for i := 0; i < ngen; i++ {
		nanctraj = append(nanctraj, rnanctraj[ngen-1-i])
	}

	return nanctraj
}
