package vcf

import (
	"errors"
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

func TestFeaturePositions(t *testing.T) {
	tests := []struct {
		Name   string
		Input  Feature
		Output FeaturePos
		Error  error
	}{{
		Name: "Simple",
		Input: Feature{
			Chrom:     "20",
			Pos:       14370,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		Output: FeaturePos{
			EndZero:   14369,
			StartZero: 14369,
			EndOne:    14370,
			StartOne:  14370,
		},
		Error: io.EOF,
	}, {
		Name: "Zero",
		Input: Feature{
			Chrom:     "20",
			Pos:       0,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		Output: FeaturePos{
			EndZero:   0,
			StartZero: 0,
			EndOne:    0,
			StartOne:  0,
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

func TestFeature_SingleGenotype(t *testing.T) {
	tests := []struct {
		Name    string
		Input   Feature
		GTOrder map[string]int
		Output  Genotype
		Error   error
	}{{
		Name: "ValidGenotype",
		Input: Feature{
			Chrom:     "20",
			Pos:       14370,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		GTOrder: map[string]int{"NA0001": 0},
		Output: Genotype{
			Id:       "NA0001",
			GT:       []int{0, 0},
			PhasedGT: true,
			Fields:   map[string]string{"GT": "0|0", "GQ": "48", "DP": "1", "HQ": "51,51"},
		},
		Error: nil,
	}, {
		Name: "MissingGenotype",
		Input: Feature{
			Chrom:     "20",
			Pos:       14370,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		GTOrder: map[string]int{"NA0002": 0},
		Error:   errors.New("genotype not in vcf"),
	}, {
		Name: "TooManyFormatLines",
		Input: Feature{
			Chrom:     "20",
			Pos:       14370,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 58, 48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		GTOrder: map[string]int{"NA0001": 0},
		Error:   errors.New("genotype has improperly formatted data"),
	}, {
		Name: "AlreadyParsed",
		Input: Feature{
			Chrom:     "20",
			Pos:       14370,
			Id:        "trs6054257",
			Ref:       "G",
			Alt:       []string{"A"},
			Qual:      29,
			Filter:    "PASS",
			Info:      map[string]string{"NS": "3", "DP": "14", "AF": "0.5", "DB": "DB", "H2": "H2"},
			Format:    map[string]int{"DP": 2, "GQ": 1, "GT": 0, "HQ": 3},
			Genotypes: [][]byte{{48, 124, 48, 58, 52, 56, 58, 49, 58, 53, 49, 44, 53, 49}},
		},
		GTOrder: map[string]int{"NA0001": 0},
		Output: Genotype{
			Id:       "NA0001",
			GT:       []int{0, 0},
			PhasedGT: true,
			Fields:   map[string]string{"GT": "0|0", "GQ": "48", "DP": "1", "HQ": "51,51"},
		},
		Error: nil,
	}}
	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			r, err := tt.Input.SingleGenotype("NA0001", tt.GTOrder)
			if !reflect.DeepEqual(err, tt.Error) {
				t.Errorf("SingleGenotype() error: unexpected error\ngot \t%v\nwant \t%v", err, tt.Error)
			} else if tt.Name == "AlreadyParsed" {
				f, _ := tt.Input.SingleGenotype("NA0001", tt.GTOrder)
				if f != r {
					t.Errorf("SingleGenotype() error: wrong pointer \ngot \t%v\nwant \t%v", f, r)
				}
			} else if r != nil { // Quick and dirty way to check there is an Output, as
				out := reflect.ValueOf(*r)
				exp := reflect.ValueOf(tt.Output)
				for i := 0; i < out.NumField(); i++ {
					got := out.Field(i).Interface()
					expt := exp.Field(i).Interface()
					if !reflect.DeepEqual(got, expt) {
						t.Errorf("SingleGenotype() error: unexpected value, check FieldOrder\ngot \t%v\nwant \t%v", got, expt)
					}
				}
				if par := map[string]*Genotype{"NA0001": r}; !reflect.DeepEqual(par, tt.Input.ParsedGenotypes) {
					t.Errorf("SingleGenotype() error: ParsedGenotypes didn't update\ngot \t%v\nwant \t%v", tt.Input.ParsedGenotypes, par)
				}
			}
		})
	}
}
