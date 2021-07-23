package unicell

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"os"
)

func DOT_Genealogy(dotfilename, popfilename string, ngen, npop int) []float64 { //Dumps genealogy of population for an epoch into a dot file, going backwards in time. Returns proportion of reproducing population
	var npars int
	var id, dadid, momid string
	indiv := NewIndiv(0)

	proptraj := make([]float64,npop)
	pop := NewPopulation(npop)
	genfile := fmt.Sprintf("%s.dot",dotfilename)
	fout, err := os.OpenFile(genfile, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(fout,"digraph G {")
	err = fout.Close()
	if err != nil {
		log.Fatal(err)
	}
	kids := make(map[int]bool)
	for i:=0; i<npop; i++ {
		kids[i]=true
	}
	pars := make(map[int]bool)
	for gen := ngen; gen>0; gen-- {
		jfilename := fmt.Sprintf("%s_%d.json",popfilename,gen)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}
		byteValue, _ := ioutil.ReadAll(popin)
		err = json.Unmarshal(byteValue, &pop)
		if err != nil {
			log.Fatal(err)
		}
		err = popin.Close()
		if err != nil{
			log.Fatal(err)
		}
		
		fout, err = os.OpenFile(genfile, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0644)
		if err != nil {
			log.Fatal(err)
		}

		for i := range(kids){
			indiv = pop.Indivs[i] 
			pars[indiv.DadId]=true
			pars[indiv.MomId]=true
			id = fmt.Sprintf("g%d:id%d",pop.Gen,indiv.Id)
			dadid = fmt.Sprintf("g%d:id%d",pop.Gen-1,indiv.DadId)
			momid = fmt.Sprintf("g%d:id%d",pop.Gen-1,indiv.MomId)
			fmt.Fprintf(fout,"\t %s-> {%s, %s}\n",id,dadid,momid)
		}
		err = fout.Close()
			if err != nil {
				log.Fatal(err)
			}
		npars = len(pars)
		proptraj = append(proptraj, float64(npars)/float64(npop))
		kids = pars //update; go back in time
		pars = make(map[int]bool) //re-initialize
	}
	return proptraj
}