package gff

import (
	"io"
	"reflect"
	"testing"
)

type FeaturePos struct {
	EndZero   uint64
	StartZero uint64
	EndOne    uint64
	StartOne  uint64
}

func TestFeature(t *testing.T) {
	tests := []struct {
		Name   string
		Input  Feature
		Output FeaturePos
		Error  error
	}{{
		Name: "ChangeOfBaseTest",
		Input: Feature{
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
		Output: FeaturePos{
			EndZero:   6484,
			StartZero: 6451,
			EndOne:    6485,
			StartOne:  6452,
		},
		Error: io.EOF,
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			res := FeaturePos{
				EndZero:   tt.Input.EndZero(),
				StartZero: tt.Input.StartZero(),
				EndOne:    tt.Input.EndOne(),
				StartOne:  tt.Input.StartOne(),
			}

			if !reflect.DeepEqual(res, tt.Output) {
				t.Errorf("feature error: unexpected read\ngot \t%v\nwant \t%v", res, tt.Output)
			}
		})
	}
}
