// Feature describes a single row of a gff3 file
// http://www.sequenceontology.org/gff3.shtml

package gff

import (
	"bytes"
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
)

// A Feature that represents a single line of a gff3 file
//
// By specification, coordinates are one-based, and any undefined
// field is to be replaced with the "." character.
type Feature struct {

	// ID of the landmark used to establish the coordinate system of the feature.
	// May contain any characters, but must escape anything not in [a-zA-Z0-9.:^*$@!+_?-|]
	Seqid string

	// Free text qualifier intended to describe the algorithm
	// or operating procedure that generated this feature.
	Source string

	// The type of the feature (previously called the "method").
	// Supposed to be constrained to a term from Sequence Ontology or an
	// SO accession number, but this check is not supported.
	Type string

	// The start of feature given in positive 1-based integer coordinates
	// relative to the Seqid in column one.
	// if "." treated as 0
	Start uint64

	// The end of feature given in positive 1-based integer coordinates
	// relative to the Seqid in column one. End must be larger than Start.
	// uf "." treated as 0
	End uint64

	// A floating point number.
	// recommended to be an E-value or P-value
	// if "." treated as math.MaxFloat64
	Score float64

	// Strand relative to landmark, one of [+,-,?]
	Strand string

	// For use with "CDS" type features, one of [0,1,2].
	// if "." treated as 3
	Phase int8

	// A semicolon separated list of <tag>=<value> pairs.
	Attributes map[string]string
}

// Column s6 (score) allows for an undefined value "."
const MissingScoreField = math.MaxFloat64

// Column  8 (score) allows for an undefined value "."
const MissingPhaseField = 3

// StartZero returns Feature.Start in zero based coordinate systems
func (f *Feature) StartZero() uint64 {
	return f.Start - 1
}

// EndZero returns Feature.End in zero based coordinate systems
func (f *Feature) EndZero() uint64 {
	return f.End - 1
}

// StartOne returns Feature.Start in one based coordinate systems (gff3 spec default)
func (f *Feature) StartOne() uint64 {
	return f.Start
}

// EndOne returns Feature.End in one based coordinate systems (gff3 spec default)
func (f *Feature) EndOne() uint64 {
	return f.End
}

// String returns the string representation of the gff3 feature
func (f *Feature) String() string {
	var start, end, score, phase, attributes string
	start = strconv.FormatUint(f.Start, 10)
	end = "."

	if f.Score == MissingScoreField {
		score = "."
	} else {
		score = strconv.FormatFloat(f.Score, 'e', -1, 64)
	}

	if p := strconv.Itoa(int(f.Phase)); p == strconv.Itoa(MissingPhaseField) {
		phase = "."
	} else {
		phase = p
	}

	if len(f.Attributes) == 0 { //Attributes is an optional column
		return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", f.Seqid, f.Source, f.Type, start, end, score, f.Strand, phase)
	} else {
		b := new(bytes.Buffer)
		var k []string
		for key := range f.Attributes {
			k = append(k, key)
		}
		sort.Strings(k)
		for _, key := range k {
			_, _ = fmt.Fprintf(b, "%s=%s;", key, f.Attributes[key])
		}
		attributes = b.String()
		attributes = strings.TrimRight(attributes, ";")
		return fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", f.Seqid, f.Source, f.Type, start, end, score, f.Strand, phase, attributes)
	}
}
