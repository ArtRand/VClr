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

func (self *VcAlignment) GroupByStrand() map[string]*VcAlignment {
	grouped := make(map[string]*VcAlignment)
	for _, r := range self.Records {
		strand := r.strand
		_, contains := grouped[strand]
		if !contains {
			vca := VcAlignmentConstruct()
			vca.AddRecord(r)
			grouped[strand] = vca
		} else {
			grouped[strand].AddRecord(r)
		}
	}
	return grouped
}

func reverseComplementBase(b string) string {
	switch {
	case b == "A":
		return "T"
	case b == "C":
		return "G"
	case b == "G":
		return "C"
	case b == "T":
		return "A"
	default:
		err := fmt.Sprintf("reverseComplementBase: Not proper input %v", b)
		panic(err)
	}
}

func correctBaseForStrand(base, strand string, forward bool) string {
	var isTemplate bool = strand == "t"
	if (isTemplate && forward) || (!isTemplate && !forward) {
		return base
	} else {
		return reverseComplementBase(base)
	}
}

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

// CallSiteOnStrand does not correct for forward/backward template/complement, it just calls the base with the argmax
// probability
func (self *VcAlignment) CallSiteOnStrand(threshold float64) string {
	site := self.Records[0].refPos
	probs := make(map[string]float64)
	call := ""
	for _, r := range self.Records {
		if r.refPos != site {
			panic("CallSiteOnStrand: Not sorted by site")
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
	if len(probs) == 0 {
		return call
	}
	normalizeProbs(&probs)
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

// CallSiteOnCodingStrand respects that there can be template and complement alignments, it corrects to the forward/
// template 'coding' orientation it aggregates the probabilities from both template and complement reads (assuming they
// are above the threshold)
func (self *VcAlignment) CallSiteOnCodingStrand(threshold float64) string {
	site := self.Records[0].refPos
	probs := make(map[string]float64)
	call := ""
	for _, r := range self.Records {
		if r.refPos != site {
			panic("CallSiteOnCodingStrand: Not sorted by site")
			return ""
		}
		if r.prob >= threshold {
			base := correctBaseForStrand(r.base, r.strand, r.forward)
			probs[base] += r.prob
		} else {
			continue
		}
	}
	if len(probs) == 0 {
		return call
	}
	normalizeProbs(&probs)
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
		siteCall , check1 := calls[site]
		rcSiteCall , check2 := calls[rcSite]
		if !check1 || !check2 {
			//fmt.Fprintln(os.Stderr, "skipping ", site, rcSite)
			//fmt.Fprintln(os.Stderr, calls)
			continue
		}
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

func CallSingleMoleculeGatcMethylation(alignment *VcAlignment, threshold float64) [][]*VariantCall {
	results := make([][]*VariantCall, 0)
	// alignment is not sorted by read, so sort by read (single molecules) first
	byRead := alignment.GroupByRead()
	for readLabel, aln := range byRead {
		// get the score for this read
		readScore := aln.ScoreRead()
		// call each site
		bySite := aln.GroupBySite()
		// strandCalls is a map of sites to calls, map[site]call
		strandCalls := make(map[int]string)
		for site, alignedPairs := range bySite {
			call := alignedPairs.CallSiteOnStrand(threshold)
			strandCalls[site] = call
		}
		sK := SortedKeys(strandCalls)
		calls := CallGatcMotifs(sK, strandCalls, readLabel, readScore)
		results = append(results, calls)
	}
	return results
}

func CallSingleMoleculeCanonicalVariants(alignment *VcAlignment, threshold float64) [][]*VariantCall {
	results := make([][]*VariantCall, 0)
	// alignment is not sorted by read, so sort by read (single molecules) first
	byRead := alignment.GroupByRead()
	for readLabel, aln := range byRead {
		// get the score for this read
		readScore := aln.ScoreRead()
		// group by site
		bySite := aln.GroupBySite()
		// strandCalls is a map of sites to calls, map[site]call
		strandCalls := make(map[int]string)
		for site, alignedPairs := range bySite {
			// call the reference position
			call := alignedPairs.CallSiteOnCodingStrand(threshold)
			strandCalls[site] = call
		}
		calls := make([]*VariantCall, 0)  // could make this length known
		for site, call := range strandCalls {
			vc := VariantCallConstruct(site, call, readLabel, readScore)
			calls = append(calls, vc)
		}
		results = append(results, calls)
	}
	return results
}

func CallSingleMoleculeMethylation(alignment *VcAlignment, threshold float64) [][]*VariantCall {
	results := make([][]*VariantCall, 0)
	// alignment is not sorted by read, so sort by read (single molecules) first
	byRead := alignment.GroupByRead()
	for readLabel, aln := range byRead {
		// get the score for this read
		readScore := aln.ScoreRead()
		// group by site
		bySite := aln.GroupBySite()
		// strandCalls is a map of sites to calls, map[site]call
		strandCalls := make(map[int]string)
		for site, alignedPairs := range bySite {
			// call the reference position
			call := alignedPairs.CallSiteOnStrand(threshold)
			strandCalls[site] = call
		}
		calls := make([]*VariantCall, 0)  // could make this length known
		for site, call := range strandCalls {
			vc := VariantCallConstruct(site, call, readLabel, readScore)
			calls = append(calls, vc)
		}
		results = append(results, calls)
	}
	return results
}

func CallSiteMethylation(siteSorted *VcAlignment, threshold float64) (string, int) {
	call := siteSorted.CallSiteOnStrand(threshold)
	coverage := coverage(siteSorted)
	return call, coverage
}

func coverage(siteSorted *VcAlignment) int {
	byRead := siteSorted.GroupByRead()
	return len(byRead)
}

func CallSite(siteSorted *VcAlignment, threshold float64) (string, int) {
	call := siteSorted.CallSiteOnCodingStrand(threshold)
	coverage := coverage(siteSorted)
	return call, coverage
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

