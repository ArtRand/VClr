package main

import (
	"flag"
	"fmt"
	vclr "github.com/ArtRand/VClr/lib"
	"math"
	"os"
	"bufio"
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

func compareCallsToReference(results [][]*vclr.VariantCall, reference string) (float64, float64) {
	if len(results) > 1 {
		err := fmt.Sprintf("compareCallsToReference: got more than one read's worth of results? %v", results)
		panic(err)
	}
	readResult := results[0]
	var numCorrect float64 = 0.0
	var numCalled float64 = 0.0
	var readScore float64 = 0.0
	for _, siteCall := range readResult {
		correctBase := reference[siteCall.RefPos]
		calledBase := siteCall.Call
		readScore = siteCall.ReadScore
		if string(correctBase) == calledBase {
			numCorrect += 1
			numCalled += 1
		} else {
			numCalled += 1
		}
	}
	return numCorrect / numCalled * 100, readScore
}

func calculatePercentCalledMethyl(results [][]*vclr.VariantCall) (float64, float64) {
	if len(results) > 1 {
		err := fmt.Sprintf("calculatePercentCalledMethyl: got more than one read's worth of results? %v", results)
		panic(err)
	}
	readResult := results[0]
	var numCalledMethyl float64 = 0.0
	var numCalled float64 = 0.0
	var readScore float64 = 0.0
	for _, siteCall := range readResult {
		readScore = siteCall.ReadScore
		if siteCall.Call == "I" || siteCall.Call == "E" {
			numCalledMethyl += 1
			numCalled += 1
		} else {
			numCalled += 1
		}
	}
	return numCalledMethyl / numCalled * 100, readScore
}

func meanFloatSlice(slc *[]float64) float64 {
	tot := 0.0
	for _, i := range *slc {
		tot += i
	}
	return tot / float64(len(*slc))
}

func callSingleStrandVariants(vca *vclr.VcAlignment, threshold *float64, reference string) {
	// first group the alignment by read
	byRead := vca.GroupByRead()
	// then get the accuracy for each strand
	templateAccuracies := make([]float64, 0)
	complementAccuracies := make([]float64, 0)
	for read, aln := range byRead {
		byStrand := aln.GroupByStrand()
		_, hasTemplate := byStrand["t"]
		_, hasComplement := byStrand["c"]
		templateAccuracy := math.NaN()
		complementAccuracy := math.NaN()
		temScore := math.NaN()
		comScore := math.NaN()
		if hasTemplate {
			templateResults := vclr.CallSingleMoleculeCanonicalVariants(byStrand["t"], *threshold)
			templateAccuracy, temScore = compareCallsToReference(templateResults, reference)
			templateAccuracies = append(templateAccuracies, templateAccuracy)
		}
		if hasComplement {
			complementResults := vclr.CallSingleMoleculeCanonicalVariants(byStrand["c"], *threshold)
			complementAccuracy, comScore = compareCallsToReference(complementResults, reference)
			complementAccuracies = append(complementAccuracies, complementAccuracy)
		}
		fmt.Fprintf(os.Stdout, "%v\t%v\t%v\t%v\t%v\n", read, templateAccuracy, complementAccuracy, temScore, comScore)
	}
	fmt.Fprintf(os.Stderr, "mean template accuracy %v\n", meanFloatSlice(&templateAccuracies))
	fmt.Fprintf(os.Stderr, "mean complement accuracy %v\n", meanFloatSlice(&complementAccuracies))
}

func callSingleStrandMethylation(vca *vclr.VcAlignment, threshold *float64) {
	// first group the alignment by read
	byRead := vca.GroupByRead()
	tempalteMethylPercents := make([]float64, 0)
	complementMethylPercents := make([]float64, 0)
	for read, aln := range byRead {
		byStrand := aln.GroupByStrand()
		_, hasTemplate := byStrand["t"]
		_, hasComplement := byStrand["c"]
		tem_percentMethyl := math.NaN()
		com_percentMethyl := math.NaN()
		temScore := math.NaN()
		comScore := math.NaN()
		if hasTemplate {
			templateResults := vclr.CallSingleMoleculeMethylation(byStrand["t"], *threshold)
			tem_percentMethyl, temScore = calculatePercentCalledMethyl(templateResults)
			tempalteMethylPercents = append(tempalteMethylPercents, tem_percentMethyl)
		}
		if hasComplement {
			complementResults := vclr.CallSingleMoleculeMethylation(byStrand["c"], *threshold)
			com_percentMethyl, comScore = calculatePercentCalledMethyl(complementResults)
			complementMethylPercents = append(complementMethylPercents, com_percentMethyl)
		}
		fmt.Fprintf(os.Stdout, "%v\t%v\t%v\t%v\t%v\n", read, tem_percentMethyl, com_percentMethyl, temScore, comScore)
	}
	fmt.Fprintf(os.Stderr, "template %% called methyl %v\n", meanFloatSlice(&tempalteMethylPercents))
	fmt.Fprintf(os.Stderr, "complement %% called methyl %v\n", meanFloatSlice(&complementMethylPercents))
}

func callSites(vca *vclr.VcAlignment, threshold *float64, canonical bool) {
	// group the alignment by site
	bySite := vca.GroupBySite()
	for site, aln := range bySite {
		var call string
		var coverage int
		if !canonical {
			call, coverage = vclr.CallSiteMethylation(aln, *threshold)
			fmt.Println(site, call, coverage)
		} else {
			call, coverage = vclr.CallSite(aln, *threshold)
			fmt.Println(site, call, coverage)
		}
	}
}

func check(ok error, msg string) {
	if ok != nil {
		panic(msg)
	} else {
		return
	}
}

func main() {
	inFile := flag.String("f", "", "file location")
	refFasta := flag.String("r", "", "reference location")
	threshold := flag.Float64("t", 0.0, "threshold")
	readScoreT := flag.Float64("s", 0.0, "readScore threshold")
	templateOnly := flag.Bool("templateOnly", false, "flag to only use template strands")
	tool := flag.String("tool", "smVariant", "Tool to use options are: \n" +
		" single molecule variant: sm-variant\n" +
		" single molecule methylation: sm-methyl\n" +
		" variant call: variant\n " +
		" methylation call: methyl")
	flag.Parse()
	vca := vclr.ParseAlignmentFile(*inFile)
	var alns *vclr.VcAlignment
	if *templateOnly {
		byStrand := vca.GroupByStrand()
		_, check := byStrand["t"]
		if !check {
			panic("Didn't find any template reads?")
		}
		alns = byStrand["t"]
	} else if *readScoreT > 0.0 {
		alns = vca.FilterByReadScore(*readScoreT)
	} else {
		alns = vca
	}

	if *tool == "sm-variant" {
		fH, ok := os.Open(*refFasta)
		check(ok, fmt.Sprintf("Error opening file %v", *inFile))
		defer fH.Close()
		fqr := vclr.FqReader{Reader: bufio.NewReader(fH)}
		r, _ := fqr.Iter()
		callSingleStrandVariants(alns, threshold, r.Seq)
	} else if *tool == "sm-methyl" {
		callSingleStrandMethylation(alns, threshold)
	} else {
		if *tool == "variant" {
			callSites(alns, threshold, true)
		} else if *tool == "methyl" {
			callSites(alns, threshold, false)
		} else {
			err := fmt.Sprintf("Error, tool %v not recognised", *tool)
			panic(err)
		}
	}

}
