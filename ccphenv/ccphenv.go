package main

// practice using gonum

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"log"
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

	var svd mat.SVD
	ok := svd.Factorize(a, mat.SVDFull)
	if !ok {
		log.Fatal("SVD failed.")
	}
	U := mat.NewDense(5, 5, nil)
	V := mat.NewDense(5, 5, nil)
	svd.UTo(U)
	svd.VTo(V)
	vals := svd.Values(nil)
	fmt.Println("Values:", vals)

	for i := 0; i < 5; i++ {
		for j := 0; j <= i; j++ {
			u1 := V.ColView(i)
			u2 := V.ColView(j)
			dot := mat.Dot(u1, u2)
			fmt.Printf("\t%e", dot)
		}
		fmt.Println()
	}
}
