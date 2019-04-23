package gff

import (
	"bytes"
	"math"
	"testing"
)

func TestWrite(t *testing.T) {
	tests := []struct {
		Name   string
		Input  Feature
		Output string
	}{{
		Name: "Full",
		Input: Feature{
			Seqid:      "Scaffold_102",
			Source:     "EVM",
			Type:       "CDS",
			Start:      6452,
			End:        6485,
			Score:      1.1e+20,
			Strand:     "+",
			Phase:      2,
			Attributes: map[string]string{"ID": "CDS705", "Parent": "mRNA906"},
		},
		Output: "##gff-version 3.2.1\nScaffold_102\tEVM\tCDS\t6452\t6485\t1.1e+20\t+\t2\tID=CDS705;Parent=mRNA906\n",
	}, {
		Name: "Empty",
		Input: Feature{
			Seqid:      ".",
			Source:     ".",
			Type:       ".",
			Start:      0,
			End:        0,
			Score:      math.MaxFloat64,
			Strand:     ".",
			Phase:      3,
			Attributes: map[string]string{},
		},
		Output: "##gff-version 3.2.1\n.\t.\t.\t.\t.\t.\t.\t.\n",
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			var b bytes.Buffer
			r, _ := NewWriter(&b)
			r.WriteFeature(&tt.Input)
			got := b.String()
			if got != tt.Output {
				t.Errorf("WriteFeature() error:\ngot \n%v want \n%v", got, tt.Output)
			}
		})
	}
}

func TestWriteAll(t *testing.T) {
	tests := []struct {
		Name   string
		Input  []Feature
		Output string
	}{{
		Name: "Full",
		Input: []Feature{
			{
				Seqid:      "Scaffold_102",
				Source:     "EVM",
				Type:       "CDS",
				Start:      6452,
				End:        6485,
				Score:      1e+20,
				Strand:     "+",
				Phase:      2,
				Attributes: map[string]string{"ID": "CDS705", "Parent": "mRNA906"},
			},
		},
		Output: "##gff-version 3.2.1\nScaffold_102\tEVM\tCDS\t6452\t6485\t1e+20\t+\t2\tID=CDS705;Parent=mRNA906\n",
	}, {
		Name: "Empty",
		Input: []Feature{
			{
				Seqid:      ".",
				Source:     ".",
				Type:       ".",
				Start:      0,
				End:        0,
				Score:      math.MaxFloat64,
				Strand:     ".",
				Phase:      3,
				Attributes: map[string]string{},
			},
		},
		Output: "##gff-version 3.2.1\n.\t.\t.\t.\t.\t.\t.\t.\n",
	}, {
		Name: "MultipleFields",
		Input: []Feature{
			{
				Seqid:      "Scaffold_103",
				Source:     "EVM",
				Type:       "CDS",
				Start:      6452,
				End:        6485,
				Score:      math.MaxFloat64,
				Strand:     "+",
				Phase:      2,
				Attributes: map[string]string{"ID": "CDS705.1", "Parent": "mRNA906"},
			},
			{
				Seqid:      "Scaffold_102",
				Source:     "EVM",
				Type:       "CDS",
				Start:      6452,
				End:        6485,
				Score:      math.MaxFloat64,
				Strand:     "+",
				Phase:      2,
				Attributes: map[string]string{"ID": "CDS705.2", "Parent": "mRNA906"},
			},
		},
		Output: "##gff-version 3.2.1\nScaffold_103\tEVM\tCDS\t6452\t6485\t.\t+\t2\tID=CDS705.1;Parent=mRNA906\nScaffold_102\tEVM\tCDS\t6452\t6485\t.\t+\t2\tID=CDS705.2;Parent=mRNA906\n",
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			var b bytes.Buffer
			r, _ := NewWriter(&b)
			var in []*Feature
			for i := range tt.Input {
				in = append(in, &tt.Input[i])
			}
			r.WriteAll(in)
			got := b.String()
			if got != tt.Output {
				t.Errorf("WriteFeature() error:\ngot \n%v want \n%v", got, tt.Output)
			}
		})
	}
}
