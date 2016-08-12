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

	for i := range results {
		for _, x := range i {
			switch {
			case x.Call == "methylated":
				methylCalls += 1
			case x.Call == "unmethylated":
				unmethylCalls += 1
			case x.Call == "hemi-methylated":
				hemiMethylCalls += 1
			}

		}
	}
	assert.True(t, unmethylCalls == 1)
	assert.True(t, methylCalls == 1)
	assert.True(t, hemiMethylCalls == 1)
}