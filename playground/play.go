package main

import (
	//"encoding/json"
	//"flag"
	"fmt"
	//"io/ioutil"
	//"log"
	//"os"
	"time"

	//"github.com/arkinjo/evodevo/multicell"
)

var T_Filename string = "traj"
var PG_Filename string  //Dump for phenotypes and genotypes
var Gid_Filename string //Genealogy of ID's
var nancfilename string
var json_in string //JSON encoding of initial population; default to empty string
var json_out string = "popout"

var CopyequalGs float64 //bugtesting variable
var JSONequalGs float64 //bugtesting variable
var DevequalGs float64  //bugtesting variable

//var test bool = false //false : training mode, true : testing mode
//This file is for playing around; testing behavior of various functions in Go etc.

/*
Objectives:
Fix genome projection bug
Test json population import/export with(out) flag O_Append

*/

func main() {
	t0 := time.Now()
	fmt.Println("Hello, world!")
	dt := time.Since(t0)
	fmt.Println("Total time taken : ",dt)
}
