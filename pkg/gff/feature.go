// Feature describes a single row of a gff3 file
// http://www.sequenceontology.org/gff3.shtml

package gff

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

// Explicit start in zero  base coordinate systems
func (f *Feature) StartZero() uint64 {
	return f.Start - 1
}

// Explicit end in zero  base coordinate systems
func (f *Feature) EndZero() uint64 {
	return f.End - 1
}

// Explicit start in one based coordinate systems (gff3 spec default)
func (f *Feature) StartOne() uint64 {
	return f.Start
}

// Explicit end in one based coordinate systems (gff3 spec default)
func (f *Feature) EndOne() uint64 {
	return f.End
}
