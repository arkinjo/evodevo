package main

import (
	"bufio"
	"encoding/json"
	"io/ioutil"
	"strconv"

	//"encoding/binary"
	"fmt"
	"log"
	"os"

	//	"fmt"
	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"

var json_in string = "upop20211223train"
var vfilename string = "upop20211223trainpvec"

//var pcafilename string = "pcatestphen5.dat"

func main() {
	var counter, prdir, r, cell, trait int
	var float, truecorr, trimmedcorr, acc float64
	var str string

	pop0 := multicell.NewPopulation(multicell.GetNcells(), multicell.GetMaxPop()) //with randomized genome to start
	multicell.SetNcells(1)

	fmt.Println("Hello, world!")
	floats := make([]float64, 0)

	file, err := os.Open("upop20211223trainpvec.dat")
	if err != nil {
		log.Fatal(err)
	}
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		counter++
		str = scanner.Text()
		if float, err = strconv.ParseFloat(str, 64); err == nil {
			//fmt.Printf("Entry:%d\tValue:%f\n", counter, float)
			//fmt.Printf("Entry:%T\tValue:%T\n", counter, float)
			floats = append(floats, float)
		}
	}
	pcacues, pcavecs := multicell.PCAtoCue(vfilename)
	for i, t := range floats {
		prdir = i / (multicell.GetNenv() * multicell.GetNcells()) //take advantage of integer division
		r = i % (multicell.GetNcells() * multicell.GetNenv())     //remainder
		cell = r / multicell.GetNenv()
		trait = r % multicell.GetNenv()
		fmt.Printf("Entry:(%d,%d,%d) , Input:%f , Output:%f\n", prdir, cell, trait, t, pcacues[prdir][cell][trait])
	}

	if json_in != "" { //read input population as a json file, if given
		fmt.Println("Importing initial population")
		jfilename := fmt.Sprintf("../pops/%s.json", json_in)
		popin, err := os.Open(jfilename)
		if err != nil {
			log.Fatal(err)
		}

		byteValue, _ := ioutil.ReadAll(popin)
		pop0.ClearGenome() //Clear genome before unmarshalling json
		err = json.Unmarshal(byteValue, &pop0)
		if err != nil {
			log.Fatal(err)
		}

		err = popin.Close()
		if err != nil {
			log.Fatal(err)
		}
		fmt.Println("Successfully imported population")
	}
	AncEnvs := pop0.RefEnvs

	for i, vec := range pcacues {
		//fmt.Printf("PCA Index: %d\t VecLength: %d \n", i, len(vec))
		truecorr = 0
		trimmedcorr = 0
		acc = 0
		for j, cell := range vec {
			//fmt.Println("Anc env :", AncEnvs[j])
			//fmt.Println("PCA env :", cell)
			//fmt.Println("PCA vec :", pcavecs[i][j])
			trimmedcorr += multicell.Innerproduct(multicell.GetTrait(AncEnvs[j]), multicell.GetTrait(cell))
			truecorr += multicell.Innerproduct(multicell.GetTrait(AncEnvs[j]), pcavecs[i][j])
			acc += multicell.DistVecs(multicell.GetTrait(cell), pcavecs[i][j])
		}
		fmt.Printf("PCA Index:%d\t True Corr: %f\t Trimmed Corr: %f\t Cue_Acc: %f \n", i, truecorr, trimmedcorr, acc)
	}

	err = file.Close()
	if err != nil {
		log.Fatal(err)
	}

}
