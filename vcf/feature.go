// Data structures to describe the parts of a VCF file
// https://samtools.github.io/hts-specs/VCFv4.3.pdfhttps://samtools.github.io/hts-specs/VCFv4.3.pdf

package vcf

import (
	"bytes"
	"errors"
	"fmt"
	"strconv"
	"strings"
)

// Meta-information line object, used for structured meta fields
type Meta struct {
	FieldType   string
	Id          string
	Number      string
	Type        string
	Description string
	Url         string
	Optional    map[string]string
	FieldOrder  []string
}

// Meta-information that's a single key-value pair
// FieldType is used for key, Id is used for value
type SingleValMeta Meta

// Header of a vcf file, including meta-information
type Header struct {
	Metas      []*Meta
	Infos      []*Meta
	Filters    []*Meta
	Formats    []*Meta
	Alts       []*Meta
	Assemblies []*Meta
	Contigs    []*Meta
	Samples    []*Meta
	Pedigrees  []*Meta
	Others     []*Meta
	SingleVals []*SingleValMeta
	PrintOrder []*Meta
	FileFormat string
	Genotypes  map[string]uint64
}

// Feature represents a single line of a vcf file.
type Feature struct {
	Chrom           string
	Pos             uint64
	Id              string
	Ref             string
	Alt             []string
	Qual            float64
	QualFormat      byte
	Filter          string
	Info            map[string]string
	InfoOrder       map[string]int
	Format          map[string]int
	Genotypes       [][]byte
	ParsedGenotypes map[string]*Genotype
}

// Genotype represents a single genotype variant in a Feature
type Genotype struct {
	Id       string
	GT       []int
	PhasedGT bool
	Fields   map[string]string
}

// OptionalToString returns string representation of any meta directive field that isn't
// ID, Number, Type, Description, or URL
func (m *Meta) optionalToString() string {
	var b = bytes.Buffer{}

	fo := m.FieldOrder

	for _, opt := range fo {
		if val, ok := m.Optional[opt]; ok {
			b.WriteString(fmt.Sprintf("%s=%s,", opt, val))
		}
	}

	return b.String()
}

func (m *Meta) raw() *Meta {
	return m
}

// String returns string representation of a ##META=<ID=VALUE,...> meta directive
func (m *Meta) String() string {
	var b = bytes.Buffer{}
	if len(m.FieldOrder) == 0 { //Generate DefaultFieldOrder if there is none
		m.FieldOrder = append(m.FieldOrder, "ID", "Number", "Type", "Description", "URL")
		for _, field := range m.Optional {
			m.FieldOrder = append(m.FieldOrder, field)
		}
	}
	b.WriteString(fmt.Sprintf("##%s=<", m.FieldType))
	for _, field := range m.FieldOrder {
		switch field {
		case "ID":
			if m.Id != "" {
				b.WriteString(fmt.Sprintf("%s=%s,", field, m.Id))
			}
		case "Number":
			if m.Number != "" {
				b.WriteString(fmt.Sprintf("%s=%s,", field, m.Number))
			}
		case "Type":
			if m.Type != "" {
				b.WriteString(fmt.Sprintf("%s=%s,", field, m.Type))
			}
		case "Description":
			if m.Description != "" {
				b.WriteString(fmt.Sprintf("%s=%s,", field, m.Description))
			}
		case "URL":
			if m.Url != "" {
				b.WriteString(fmt.Sprintf("URL=%s,", m.Url))
			}
		default:
			continue
		}
	}

	if len(m.Optional) != 0 {
		b.WriteString(m.optionalToString())
	}
	str := b.String()
	str = strings.TrimRight(str, ",")
	return fmt.Sprintf("%s>", str)
}

// String returns string representation of a ##META=VALUE style meta directive
func (m *SingleValMeta) String() string {
	return fmt.Sprintf("##%s=%s", m.FieldType, m.Id)
}

func NewHeader() *Header {
	var h Header
	h.Metas = make([]*Meta, 0)
	h.Infos = make([]*Meta, 0)
	h.Filters = make([]*Meta, 0)
	h.Formats = make([]*Meta, 0)
	h.Alts = make([]*Meta, 0)
	h.Assemblies = make([]*Meta, 0)
	h.Contigs = make([]*Meta, 0)
	h.Samples = make([]*Meta, 0)
	h.Pedigrees = make([]*Meta, 0)
	h.Others = make([]*Meta, 0)
	h.SingleVals = make([]*SingleValMeta, 0)
	h.PrintOrder = make([]*Meta, 0)
	h.Genotypes = make(map[string]uint64)
	return &h
}

// StartZero returns Feature.Start in zero based coordinate systems
func (f *Feature) StartZero() uint64 {
	if f.Pos == 0 {
		return 0
	}
	return f.Pos - 1
}

// EndZero returns Feature.End in zero based coordinate systems
func (f *Feature) EndZero() uint64 {
	return f.StartZero()
}

// StartOne returns Feature.Start in one based coordinate systems (gff3 spec default)
func (f *Feature) StartOne() uint64 {
	return f.Pos
}

// EndOne returns Feature.End in one based coordinate systems (gff3 spec default)
func (f *Feature) EndOne() uint64 {
	return f.StartOne()
}

//SingleGenotype returns a pointer to a Genotype or an error
func (f *Feature) SingleGenotype(gen string, order map[string]uint64) (*Genotype, error) {
	if loc, ok := order[gen]; ok { //gen is a valid genotype
		if preParsed, ok := f.ParsedGenotypes[gen]; ok { //gen has already been accessed for this feature
			return preParsed, nil
		} else { //gen needs to be extracted from the info field
			info := bytes.Split(f.Genotypes[loc], []byte{':'})
			if len(info) != len(f.Format) { //info is improperly formatted
				return nil, errors.New("genotype has improperly formatted data")
			} else {
				parsedGT := Genotype{PhasedGT: false}
				parsedGT.Id = gen
				parsedGT.Fields = make(map[string]string, len(f.Format))
				for key, value := range f.Format {
					parsedGT.Fields[key] = string(info[value])
					if key == "GT" {
						gt := bytes.Split(info[value], []byte{'|'})
						if len(gt) > 1 {
							parsedGT.PhasedGT = true
						} else if len(gt) == 1 {
							gt = bytes.Split(info[value], []byte{'/'})
						}

						parsedGT.GT = make([]int, len(gt))
						for i := range gt {
							if string(gt[i]) == "." {
								parsedGT.GT[i] = -1
							} else {
								val, _ := strconv.Atoi(string(gt[i]))
								parsedGT.GT[i] = val
							}
						}
					}
				}
				if len(f.ParsedGenotypes) == 0 {
					f.ParsedGenotypes = make(map[string]*Genotype)
				}
				f.ParsedGenotypes[gen] = &parsedGT
				return &parsedGT, nil
			}
		}
	} else {
		return nil, errors.New("genotype not in vcf")
	}
}

//MultipleGenotypes returns an array of pointers to genotypes, along with an array of errors
func (f *Feature) MultipleGenotypes(gens []string, order map[string]uint64) ([]*Genotype, []error) {
	gts := make([]*Genotype, len(gens))
	errs := make([]error, len(gens))
	for i, gen := range gens {
		gt, err := f.SingleGenotype(gen, order)
		gts[i] = gt
		errs[i] = err
	}
	return gts, errs
}

//AllGenotypes returns an array of pointers to all genotypes, along with any errors
func (f *Feature) AllGenotypes(order map[string]uint64) ([]*Genotype, []error) {
	gts := make([]string, len(order))
	for gt, i := range order {
		gts[i] = gt
	}
	return f.MultipleGenotypes(gts, order)
}
