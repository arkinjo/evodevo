package main

import (
	"bufio"
	"strconv"

	//"encoding/binary"
	"fmt"
	"log"
	"os"

	//	"fmt"
	"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"

//var jfilename string = "jsonpcatest.json"
var vfilename string = "m2pcavectest.dat"

//var pcafilename string = "pcatestphen5.dat"

func main() {
	var counter, prdir, r, cell, trait int
	var float float64
	var str string

	multicell.SetNcells(2)

	fmt.Println("Hello, world!")
	floats := make([]float64, 0)

	file, err := os.Open(vfilename)
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
	pcacues := multicell.PCAtoCue(vfilename)
	for i, t := range floats {
		prdir = i / (multicell.GetNenv() * multicell.GetNcells()) //take advantage of integer division
		r = i % (multicell.GetNcells() * multicell.GetNenv())     //remainder
		cell = r / multicell.GetNenv()
		trait = r % multicell.GetNenv()
		fmt.Printf("Entry:(%d,%d,%d) , Input:%f , Output:%f\n", prdir, cell, trait, t, pcacues[prdir][cell][trait])
	}

	err = file.Close()
	if err != nil {
		log.Fatal(err)
	}

}
