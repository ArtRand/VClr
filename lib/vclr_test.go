package VClr

import (
	"testing"
	"github.com/stretchr/testify/assert"
	"os"
	"fmt"
)

func getTestingFile() string {
	file := "../oneOfEach.tsv"
	_, ok := os.Stat(file)
	if ok != nil {
		panic("Didn't find testing file")
	}
	return file
}

func TestParseAlignmentFile(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	assert.True(t, len(vca.Records) == 70, "Incorrect number of records, got %v", len(vca.Records))
}

func TestVcAlignment_GroupByRead(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	byRead := vca.GroupByRead()
	reads := 0
	for _ = range byRead {
		reads += 1
	}
	assert.True(t, len(byRead) == 3)
}

func TestVcAlignment_GroupBySite(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	bySite := vca.GroupBySite()
	assert.True(t, len(bySite) == 2)
}

func TestVcAlignment_GroupByStrand(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	byStrand := vca.GroupByStrand()
	for strand, aln := range byStrand {
		for _, r := range aln.Records {
			assert.True(t, r.strand == strand)
		}
	}
}

func TestCallGatcMotifs(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	results := CallSingleMoleculeGatcMethylation(vca, 0.1)
	unmethylCalls := 0
	methylCalls := 0
	hemiMethylCalls := 0

	for _, read := range results {
		for _, site := range read {
			switch {
			case site.Call == "methylated":
				methylCalls += 1
			case site.Call == "unmethylated":
				unmethylCalls += 1
			case site.Call == "hemi-methylated":
				hemiMethylCalls += 1
			}

		}
	}
	assert.True(t, unmethylCalls == 1)
	assert.True(t, methylCalls == 1)
	assert.True(t, hemiMethylCalls == 1)
}

func TestCallSingleMoleculeMethylation(t *testing.T) {
	file := getTestingFile()
	vca := ParseAlignmentFile(file)
	results := CallSingleMoleculeMethylation(vca, 0.0)
	fmt.Println(results)
}

func TestCallSingleMoleculeCanonicalVariants(t *testing.T) {
	file := "../canonical_vc.tsv"
	_, ok := os.Stat(file)
	if ok != nil {
		panic("Didn't find testing file")
	}
	vca := ParseAlignmentFile(file)
	results := CallSingleMoleculeCanonicalVariants(vca, 0.1)
	assert.True(t, len(results) == 1)
	for _, readResult := range results {
		assert.True(t, len(readResult) == 18)
		for _, siteCall := range readResult {
			assert.True(t, siteCall.Call == "A")
		}
	}
}

func strandAccuracyTest(results [][]*VariantCall) float64 {
	var correct float64 = 0.0
	var totCalls float64 = 0.0
	for _, readResult := range results {
		for _, siteCall := range readResult {
			if siteCall.Call == "A" {
				correct += 1
				totCalls += 1
			} else {
				totCalls += 1
			}
		}
	}
	return correct / totCalls * 100
}

func TestStrandAccuracy(t *testing.T) {
	file := "../canonical_vc.tsv"
	_, ok := os.Stat(file)
	if ok != nil {
		panic("Didn't find testing file")
	}
	vca := ParseAlignmentFile(file)
	byStrand := vca.GroupByStrand()
	results := CallSingleMoleculeCanonicalVariants(byStrand["t"], 0.1)
	percentCorrect := strandAccuracyTest(results)
	assert.True(t, percentCorrect >= 80)
	results = CallSingleMoleculeCanonicalVariants(byStrand["c"], 0.1)
	percentCorrect = strandAccuracyTest(results)
	assert.True(t, percentCorrect >= 80)
}

func TestCallSiteMethylation (t *testing.T) {
	file := "../probs.tsv"
	_, ok := os.Stat(file)
	if ok != nil {
		panic("Didn't find testing file")
	}
	vca := ParseAlignmentFile(file)
	bySite_template := vca.GroupByStrand()["t"].GroupBySite()
	//bySite := vca.GroupBySite()
	for site, aln := range bySite_template {
		call := CallSiteMethylation(aln, 0.0)
		fmt.Println(site, call)
	}
}