package VClr

import (
	"os"
	"fmt"
	"encoding/csv"
	"io"
	"strconv"
	"math"
	"sort"
)

type AlnRecord struct {
	refPos int
	base string
	prob float64
	strand string
	forward bool
	readLabel string
}

func (self AlnRecord) String() string {
	return fmt.Sprintf("pos:%v base:%v prob:%v strand:%v forward: %v read:%v",
		self.refPos, self.base, self.prob, self.strand, self.forward, self.readLabel)
}

func AlnRecordConstruct(refPos int, base, strand, readLabel, forwardStr string, prob float64) *AlnRecord {
	var forward bool
	if forwardStr == "forward" {
		forward = true
	} else {
		forward = false
	}
	return &AlnRecord{refPos: refPos, base: base, prob: prob, strand: strand, forward: forward, readLabel: readLabel}
}

type VcAlignment struct {
	Records []*AlnRecord
}

func VcAlignmentConstruct() *VcAlignment {
	return &VcAlignment{Records: make([]*AlnRecord, 0)}
}

func (self VcAlignment) String() string {
	s := fmt.Sprintf("%v\n", self.Records[0])
	for i := 1; i < len(self.Records); i++ {
		s += fmt.Sprintf("%v\n", self.Records[i])
	}
	return s
}

func (self *VcAlignment) AddRecord(r *AlnRecord) {
	self.Records = append(self.Records, r)
}

func (self *VcAlignment) GroupByRead() map[string]*VcAlignment {
	grouped := make(map[string]*VcAlignment)
	for _, r := range self.Records {
		rL := r.readLabel
		_, contains := grouped[rL]
		if !contains {
			vca := VcAlignmentConstruct()
			vca.AddRecord(r)
			grouped[rL] = vca
		} else {
			grouped[rL].AddRecord(r)
		}
	}
	return grouped
}

func (self *VcAlignment) GroupBySite() map[int]*VcAlignment {
	grouped := make(map[int]*VcAlignment)
	for _, r := range self.Records {
		site := r.refPos
		_, contains := grouped[site]
		if !contains {
			vca := VcAlignmentConstruct()
			vca.AddRecord(r)
			grouped[site] = vca
		} else {
			grouped[site].AddRecord(r)
		}
	}
	return grouped
}

//func reverseComplement

//func correctBaseForStrand(base, strand string, forward bool) string {
//	if
//}

func normalizeProbs(probs *map[string]float64) {
	var total float64 = 0.0
	for _, v := range *probs {
		total += v
	}
	for k := range *probs {
		(*probs)[k] /= total
	}
}

func checkProbs(probs map[string]float64, allowedDiff float64) bool {
	var total float64 = 0.0
	for _, v := range probs {
		total += v
	}
	if math.Abs(1.0 - total) > allowedDiff {
		return false
	} else {
		return true
	}
}

func (self *VcAlignment) CallSite(threshold float64) string {
	site := self.Records[0].refPos
	probs := make(map[string]float64)
	call := ""
	for _, r := range self.Records {
		if r.refPos != site {
			fmt.Println("CallSite: Not sorted by site")
			return ""
		}
		if r.prob >= threshold {
			base := r.base
			prob := r.prob
			probs[base] += prob
		} else {
			continue
		}
	}
	normalizeProbs(&probs)
	if len(probs) == 0 {
		return call
	}
	probsCheck := checkProbs(probs, 0.01)
	if !probsCheck {
		err := fmt.Sprintf("normalization didn't work probs: %v", probs)
		panic(err)
	}

	maxProb := math.Inf(-1)

	for base, prob := range probs {
		if prob > maxProb {
			maxProb = prob
			call = base
		}
	}
	return call
}

func (self *VcAlignment) ScoreRead() float64 {
	firstReadLabel := self.Records[0].readLabel
	var total float64 = 0.0
	for _, r := range self.Records {
		if r.readLabel != firstReadLabel {
			panic("ScoreRead: Not sorted by read")
		}
		total += r.prob
	}
	return 100 * total / float64(len(self.Records))
}

func SortedKeys(m map[int]string) []int {
	sK := make([]int, 0)
	for k := range m {
		sK = append(sK, k)
	}
	sort.Ints(sK)
	return sK
}

type VariantCall struct {
	RefPos int
	ReadLabel string
	ReadScore float64
	Call   string
}

func VariantCallConstruct(refPos int, call string, readLabel string, readScore float64) *VariantCall {
	return &VariantCall{RefPos: refPos, Call: call, ReadLabel: readLabel, ReadScore: readScore}
}

func (self VariantCall) String() string {
	return fmt.Sprintf("(%v - %v)", self.RefPos, self.Call)
}

// calls the motifs on each read
func CallGatcMotifs(sortedSites []int, calls map[int]string, readLabel string, readScore float64) []*VariantCall {
	variantCalls := make([]*VariantCall, 0)
	for i := 0; i < len(sortedSites); i += 2 {
		// check that we have both sites in the prob table
		site := sortedSites[i]
		rcSite := site + 1
		_, check1 := calls[site]
		_, check2 := calls[rcSite]
		if !check1 || !check2 {
			fmt.Fprintln(os.Stderr, "skipping ", site, rcSite)
			fmt.Fprintln(os.Stderr, calls)
			continue
		}
		siteCall := calls[site]
		rcSiteCall := calls[rcSite]
		switch {
		case siteCall == "" || rcSiteCall == "":
			vc := VariantCallConstruct(site, "unclassified", readLabel, readScore)
			variantCalls = append(variantCalls, vc)
		case siteCall == rcSiteCall && siteCall == "A":
			vc := VariantCallConstruct(site, "unmethylated", readLabel, readScore)
			variantCalls = append(variantCalls, vc)
		case siteCall == rcSiteCall && siteCall == "I":
			vc := VariantCallConstruct(site, "methylated", readLabel, readScore)
			variantCalls = append(variantCalls, vc)
		case siteCall != rcSiteCall:
			vc := VariantCallConstruct(site, "hemi-methylated", readLabel, readScore)
			variantCalls = append(variantCalls, vc)
		}
	}
	return variantCalls
}

func CallSingleMoleculeGatcMethylation(alignment *VcAlignment, threshold float64) <- chan []*VariantCall {
	results := make(chan []*VariantCall)
	// alignment is not sorted by read, so sort by read (single molecules) first
	byRead := alignment.GroupByRead()
	go func() {
		for readLabel, aln := range byRead {
			//go callGatcRead(aln, threshold, readLabel, results)
			// get the score for this read
			readScore := aln.ScoreRead()
			// call each site
			bySite := aln.GroupBySite()
			// strandCalls is a map of sites to calls, map[site]call
			strandCalls := make(map[int]string)
			for site, alignedPairs := range bySite {
				call := alignedPairs.CallSite(threshold)
				strandCalls[site] = call
			}
			sK := SortedKeys(strandCalls)
			calls := CallGatcMotifs(sK, strandCalls, readLabel, readScore)
			results <- calls
		}
		close(results)
	}()
	//for readLabel, aln := range byRead {
		//go callGatcRead(aln, threshold, readLabel, results)
		// get the score for this read
		//readScore := aln.ScoreRead()
		// call each site
		//bySite := aln.GroupBySite()
		// strandCalls is a map of sites to calls, map[site]call
		//strandCalls := make(map[int]string)
		//for site, alignedPairs := range bySite {
		//	call := alignedPairs.CallSite(threshold)
		//	strandCalls[site] = call
		//}
		//sK := SortedKeys(strandCalls)
		//calls := CallGatcMotifs(sK, strandCalls, readLabel, readScore)
		//results <- calls
	//}
	return results
}

func ParseAlignmentFile(filePath string) *VcAlignment {
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println(err)
		return nil
	}
	defer file.Close()
	r := csv.NewReader(file)
	r.Comma = '\t'
	VcA := VcAlignmentConstruct()
	for {
		record, err := r.Read()
		// end-of-file is fitted into err
		if err == io.EOF {
			break
		} else if err != nil {
			fmt.Println("Error:", err)
			return nil
		}
		// record is an array of string so is directly printable
		refPos, _ := strconv.Atoi(record[1])
		base := record[2]
		prob, _ := strconv.ParseFloat(record[3], 64)
		strand := record[4]
		forwardStr := record[5]
		readLabel := record[6]
		aR := AlnRecordConstruct(refPos, base, strand, readLabel, forwardStr, prob)
		VcA.Records = append(VcA.Records, aR)
	}
	return VcA
}

