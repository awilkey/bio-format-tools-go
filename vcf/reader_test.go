package vcf

import (
	"errors"
	"reflect"
	"strings"
	"testing"
)

func TestNewReader(t *testing.T) {
	tests := []struct {
		Name   string
		Input  string
		Output Header
		Error  error
	}{{
		Name:  "Minimum",
		Input: "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		Output: Header{
			Metas:      make([]*Meta, 0),
			Infos:      make([]*Meta, 0),
			Filters:    make([]*Meta, 0),
			Formats:    make([]*Meta, 0),
			Alts:       make([]*Meta, 0),
			Assemblies: make([]*Meta, 0),
			Contigs:    make([]*Meta, 0),
			Samples:    make([]*Meta, 0),
			Pedigrees:  make([]*Meta, 0),
			Others:     make([]*Meta, 0),
			PrintOrder: make([]*Meta, 0),
			SingleVals: make([]*SingleValMeta, 0),
			FileFormat: "VCFv4.2",
			Genotypes:  make(map[string]uint64),
		},
	}, {
		Name:  "SingleValueMeta",
		Input: "##fileformat=VCFv4.2\n##source=myImputationProgramV3.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		Output: Header{
			Metas:      make([]*Meta, 0),
			Infos:      make([]*Meta, 0),
			Filters:    make([]*Meta, 0),
			Formats:    make([]*Meta, 0),
			Alts:       make([]*Meta, 0),
			Assemblies: make([]*Meta, 0),
			Contigs:    make([]*Meta, 0),
			Samples:    make([]*Meta, 0),
			Pedigrees:  make([]*Meta, 0),
			Others:     make([]*Meta, 0),
			SingleVals: []*SingleValMeta{
				{
					FieldType: "source",
					Id:        "myImputationProgramV3.1",
				},
			},
			PrintOrder: make([]*Meta, 0),
			FileFormat: "VCFv4.2",
			Genotypes:  make(map[string]uint64),
		},
	}, {
		Name:  "SingleComplexMeta",
		Input: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		Output: Header{
			Metas: make([]*Meta, 0),
			Infos: []*Meta{
				{
					FieldType:  "INFO",
					Id:         "myImputationProgramV3.1",
					FieldOrder: []string{"ID"},
				},
			},
			Filters:    make([]*Meta, 0),
			Formats:    make([]*Meta, 0),
			Alts:       make([]*Meta, 0),
			Assemblies: make([]*Meta, 0),
			Contigs:    make([]*Meta, 0),
			Samples:    make([]*Meta, 0),
			Pedigrees:  make([]*Meta, 0),
			Others:     make([]*Meta, 0),
			SingleVals: make([]*SingleValMeta, 0),
			PrintOrder: []*Meta{
				{
					FieldType:  "INFO",
					Id:         "myImputationProgramV3.1",
					FieldOrder: []string{"ID"},
				},
			},
			FileFormat: "VCFv4.2",
			Genotypes:  make(map[string]uint64),
		},
	}, {
		Name:  "MetaAndGenotype",
		Input: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001",
		Output: Header{
			Metas: make([]*Meta, 0),
			Infos: []*Meta{
				{
					FieldType:  "INFO",
					Id:         "myImputationProgramV3.1",
					FieldOrder: []string{"ID"},
				},
			},
			Filters:    make([]*Meta, 0),
			Formats:    make([]*Meta, 0),
			Alts:       make([]*Meta, 0),
			Assemblies: make([]*Meta, 0),
			Contigs:    make([]*Meta, 0),
			Samples:    make([]*Meta, 0),
			Pedigrees:  make([]*Meta, 0),
			Others:     make([]*Meta, 0),
			SingleVals: make([]*SingleValMeta, 0),
			PrintOrder: []*Meta{
				{
					FieldType:  "INFO",
					Id:         "myImputationProgramV3.1",
					FieldOrder: []string{"ID"},
				},
			},
			FileFormat: "VCFv4.2",
			Genotypes:  map[string]uint64{"NA00001": 0},
		},
	}, {
		Name:  "MetaAndGenotypeNoFormat",
		Input: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNA00001",
		Error: errors.New("header needs a FORMAT column before adding genotypes"),
	}, {
		Name:  "MetaAndGenotypeNoGenotype",
		Input: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
		Error: errors.New("header FORMAT column must be followed by at least one genotype"),
	}, {
		Name:  "ShortHeader",
		Input: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER",
		Error: errors.New("header has too few columns to be minimum vcf"),
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			r, err := NewReader(strings.NewReader(tt.Input))
			if !reflect.DeepEqual(err, tt.Error) {
				t.Errorf("Read() error: unexpected error\ngot \t%v\nwant \t%v", err, tt.Error)
			} else if len(tt.Output.FileFormat) > 0 { // Quick and dirty way to check there is an Output, as
				out := reflect.ValueOf(*r.Header)
				exp := reflect.ValueOf(tt.Output)
				for i := 0; i < out.NumField(); i++ {
					got := out.Field(i).Interface()
					expt := exp.Field(i).Interface()
					if !reflect.DeepEqual(got, expt) {
						t.Errorf("Read() error: unexpected value, check FieldOrder\ngot \t%v\nwant \t%v", got, expt)
					}
				}
			}
		})
	}
}
