package vcf

import (
	"bytes"
	"testing"
)

func TestWriteHeader(t *testing.T) {
	tests := []struct {
		Name   string
		Input  Header
		Output string
	}{{
		Name: "Minimum",
		Input: Header{
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
		Output: "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}, {
		Name: "SingleValueMeta",
		Input: Header{
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
		Output: "##fileformat=VCFv4.2\n##source=myImputationProgramV3.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}, {
		Name: "SingleComplexMeta",
		Input: Header{
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
			SingleVals: []*SingleValMeta{
				{
					FieldType: "source",
					Id:        "myImputationProgramV3.1",
				},
			},
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
		Output: "##fileformat=VCFv4.2\n##source=myImputationProgramV3.1\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}, {
		Name: "SingleAndMeta",
		Input: Header{
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
		Output: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
	}, {
		Name:   "MetaAndGenotype",
		Output: "##fileformat=VCFv4.2\n##INFO=<ID=myImputationProgramV3.1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001",
		Input: Header{
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
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			var b bytes.Buffer
			r, _ := NewWriter(&b)
			r.WriteHeader(tt.Input)
			got := b.String()
			if got != tt.Output {
				t.Errorf("WriteHeader() error:\ngot \n%v \nwant \n%v", got, tt.Output)
			}
		})
	}
}

func TestWriteFeature(t *testing.T) {
	tests := []struct {
		Name   string
		Input  Feature
		Output string
	}{{
		Name: "Minimum",
		Input: Feature{
			Chrom:      "20",
			Pos:        14370,
			Id:         "trs6054257",
			Ref:        "G",
			Alt:        []string{"A"},
			Qual:       29,
			QualFormat: 'f',
			Filter:     "PASS",
			Info:       map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			InfoOrder: map[string]int{
				"NS": 0,
				"DP": 1,
				"AF": 2,
				"DB": 3,
				"H2": 4,
			},
		},
		Output: "\n20\t14370\ttrs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2",
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			var b bytes.Buffer
			r, _ := NewWriter(&b)
			r.WriteFeature(&tt.Input)
			got := b.String()
			if got != tt.Output {
				t.Errorf("WriteHeader() error:\ngot \n%v \nwant \n%v", got, tt.Output)
			}
		})
	}
}
