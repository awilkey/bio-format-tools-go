package gff_test

import (
	"bio-format-tools-go/pkg/gff"
	"io"
	"math"
	"reflect"
	"strings"
	"testing"
)

func TestRead(t *testing.T) {
	tests := []struct {
		Name   string
		Input  string
		Output gff.Feature
		Error  error
	}{{
		Name: "Full",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2	ID=CDS705;Parent=mRNA906",
		Output: gff.Feature{
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
		Error: io.EOF,
	}, {
		Name: "Empty",
		Input: ".	.	.	.	.	.	.	.	.",
		Output: gff.Feature{
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
		Error: io.EOF,
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			r := gff.NewReader(strings.NewReader(tt.Input))
			out, err := r.Read()
			if !reflect.DeepEqual(err, tt.Error) {
				t.Errorf("Read() error:\ngot \t%v\nwant \t%v", err, tt.Error)
			} else if !reflect.DeepEqual(*out, tt.Output) {
				println(len(out.Attributes), len(tt.Output.Attributes))
				t.Errorf("Read() error:\ngot \t%v\nwant \t%v", *out, tt.Output)
			}
		})
	}
}
