package main

import (
	"flag"
	"fmt"
	vclr "github.com/ArtRand/VClr/VClr"
)

func main() {
	inFile := flag.String("f", "", "file location")
	threshold := flag.Float64("t", 0.0, "threshold")
	flag.Parse()
	vca := vclr.ParseAlignmentFile(*inFile)
	//fmt.Println(*threshold)
	results := vclr.CallSingleMoleculeGatcMethylation(vca, *threshold)
	for i := range results {
		var methylCalls float64 = 0.0
		var unmethylCalls float64 = 0.0
		var hemiMethylCalls float64 = 0.0
		var unclassifiedCalls float64 = 0.0
		for _, x := range i {
			switch {
			case x.Call == "methylated":
				fmt.Println(x.Call, x.ReadLabel, x.RefPos, x.ReadScore)
				methylCalls += 1
			case x.Call == "unmethylated":
				fmt.Println(x.Call, x.ReadLabel, x.RefPos, x.ReadScore)
				unmethylCalls += 1
			case x.Call == "hemi-methylated":
				fmt.Println(x.Call, x.ReadLabel, x.RefPos, x.ReadScore)
				hemiMethylCalls += 1
			case x.Call == "unclassified":
				fmt.Println(x.Call, x.ReadLabel, x.RefPos, x.ReadScore)
				unclassifiedCalls += 1
			}
		}
	}
}
