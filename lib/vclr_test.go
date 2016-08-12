package VClr

import (
	"testing"
	"github.com/stretchr/testify/assert"
	"os"
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