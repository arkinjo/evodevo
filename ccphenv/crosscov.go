package main

// practice using gonum

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math/rand"
)

func main() {
	a := mat.NewDense(5, 5, nil)
	for i := 0; i < 5; i++ {
		v := rand.NormFloat64()
		for j := 0; j < 5; j++ {
			w := rand.NormFloat64()
			a.Set(i, j, v*w)
		}
	}
	fmt.Println(mat.Trace(a))
}
