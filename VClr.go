package main

import (
	"flag"
	"fmt"
	vclr "github.com/ArtRand/VClr/lib"
)

func callGatcMethylation(vca *vclr.VcAlignment, threshold *float64) {
	results := vclr.CallSingleMoleculeGatcMethylation(vca, *threshold)
	for _, readCalls := range results {
		var methyl float64 = 0.0
		var hemi float64 = 0.0
		//var uncl float64 = 0.0
		var unmethyl float64 = 0.0
		var thisRead string = ""
		var thisScore float64
		for _, site := range readCalls {
			thisRead = site.ReadLabel
			thisScore = site.ReadScore
			switch {
			case site.Call == "methylated":
				methyl += 1
			case site.Call == "unmethylated":
				unmethyl += 1
			case site.Call == "hemi-methylated":
				hemi += 1
			}

		}
		totMethylCalls := methyl + hemi + unmethyl
		if totMethylCalls == 0 {
			//fmt.Println(thisRead, "<<Skipped?")
			continue
		}
		perMeth := 100 * methyl / totMethylCalls
		perUnmeth := 100 * unmethyl / totMethylCalls
		perHemi := 100 * hemi / totMethylCalls
		fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\n", thisRead, perUnmeth, perMeth, perHemi, totMethylCalls, thisScore)
		//fmt.Println(thisRead, perUnmeth, perMeth, perHemi, totMethylCalls)
	}
}

func main() {
	inFile := flag.String("f", "", "file location")
	threshold := flag.Float64("t", 0.0, "threshold")
	flag.Parse()
	vca := vclr.ParseAlignmentFile(*inFile)
	callGatcMethylation(vca, threshold)
}
