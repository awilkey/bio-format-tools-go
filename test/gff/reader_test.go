package gff_test

import (
	"errors"
	"github.com/awilkey/bio-format-tools-go/pkg/gff"
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
	}, {
		Name: "ShortField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2",
		Output: gff.Feature{
			Seqid:  "Scaffold_102",
			Source: "EVM",
			Type:   "CDS",
			Start:  6452,
			End:    6485,
			Score:  1e+20,
			Strand: "+",
			Phase:  2,
		},
		Error: io.EOF,
	}, {
		Name: "ErrorShortField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+",
		Error: errors.New("wrong number of fields"),
	}, {
		Name: "TooManyField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2	ID=Parent=mRNA906;	Alice=Bob;",
		Error: errors.New("wrong number of fields"),
	}, {
		Name:  "Comment",
		Input: "# This is a comment.",
		Error: io.EOF,
	}, {
		Name:  "EmptyLine",
		Input: "\n",
		Error: io.EOF,
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			r := gff.NewReader(strings.NewReader(tt.Input))
			out, err := r.Read()
			if !reflect.DeepEqual(err, tt.Error) {
				t.Errorf("Read() error: unexpected error\ngot \t%v\nwant \t%v", err, tt.Error)
			} else if out != nil && !reflect.DeepEqual(*out, tt.Output) {
				t.Errorf("Read() error: unexpected read\ngot \t%v\nwant \t%v", *out, tt.Output)
			}
		})
	}
}

func TestReadAll(t *testing.T) {
	tests := []struct {
		Name   string
		Input  string
		Output []gff.Feature
		Error  error
	}{{
		Name: "Full",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2	ID=CDS705;Parent=mRNA906",
		Output: []gff.Feature{
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
		Error: io.EOF,
	}, {
		Name: "Empty",
		Input: ".	.	.	.	.	.	.	.	.",
		Output: []gff.Feature{
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
		Error: io.EOF,
	}, {
		Name: "ShortField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2",
		Output: []gff.Feature{
			{
				Seqid:  "Scaffold_102",
				Source: "EVM",
				Type:   "CDS",
				Start:  6452,
				End:    6485,
				Score:  1e+20,
				Strand: "+",
				Phase:  2,
			},
		},
		Error: io.EOF,
	}, {
		Name: "ErrorShortField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+",
		Error: errors.New("wrong number of fields"),
	}, {
		Name: "TooManyField",
		Input: "Scaffold_102	EVM	CDS	6452	6485	1e20	+	2	ID=Parent=mRNA906;	Alice=Bob;",
		Error: errors.New("wrong number of fields"),
	}, {
		Name:  "Comment",
		Input: "# This is a comment.",
		Error: io.EOF,
	}, {
		Name:  "EmptyLine",
		Input: "\n",
		Error: io.EOF,
	}, {
		Name: "MultipleFields",
		Input: `Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.1;Parent=mRNA906
Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.2;Parent=mRNA906`,
		Output: []gff.Feature{
			{
				Seqid:      "Scaffold_102",
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
		Error: io.EOF,
	}, {
		Name: "MultipleFieldsWithComment",
		Input: `#This is a comment.
Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.1;Parent=mRNA906
Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.2;Parent=mRNA906`,
		Output: []gff.Feature{
			{
				Seqid:      "Scaffold_102",
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
		Error: io.EOF,
	}, {
		Name: "MultipleFieldsWithEmptyLine",
		Input: `Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.1;Parent=mRNA906

Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.2;Parent=mRNA906`,
		Output: []gff.Feature{
			{
				Seqid:      "Scaffold_102",
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
		Error: io.EOF,
	}, {
		Name: "CommentAndField",
		Input: `#This is a comment
Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.2;Parent=mRNA906`,
		Output: []gff.Feature{
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
		Error: io.EOF,
	}, {
		Name: "CommentAndEmptyLine",
		Input: `#This is a Comment

`,
		Error: io.EOF,
	}, {
		Name: "FullAndShortFeature",
		Input: `Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.1;Parent=mRNA906
Scaffold_102	EVM	CDS	6452	6485	.	+	2`,
		Output: []gff.Feature{
			{
				Seqid:      "Scaffold_102",
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
				Seqid:  "Scaffold_102",
				Source: "EVM",
				Type:   "CDS",
				Start:  6452,
				End:    6485,
				Score:  math.MaxFloat64,
				Strand: "+",
				Phase:  2,
			},
		},
		Error: io.EOF,
	}, {
		Name: "FullAndTooShortFeature",
		Input: `Scaffold_102	EVM	CDS	6452	6485	.	+	2	ID=CDS705.1;Parent=mRNA906
Scaffold_102	EVM	CDS	6452	6485	.	+`,
		Output: []gff.Feature{
			{
				Seqid:      "Scaffold_102",
				Source:     "EVM",
				Type:       "CDS",
				Start:      6452,
				End:        6485,
				Score:      math.MaxFloat64,
				Strand:     "+",
				Phase:      2,
				Attributes: map[string]string{"ID": "CDS705.1", "Parent": "mRNA906"},
			},
		},
		Error: errors.New("wrong number of fields"),
	}}

	for _, tt := range tests {
		t.Run(tt.Name, func(t *testing.T) {
			r := gff.NewReader(strings.NewReader(tt.Input))
			out, err := r.ReadAll()
			var res []gff.Feature
			for _, of := range out {
				res = append(res, *of)
			}

			if !reflect.DeepEqual(err, tt.Error) {
				t.Errorf("ReadAll() error: unexpected error\ngot \t%v\nwant \t%v", err, tt.Error)
			} else if len(out) != len(tt.Output) {
				t.Errorf("ReadAll() error: number of features\ngot \t%v\nwant \t%v", len(out), len(tt.Output))
			} else if !reflect.DeepEqual(res, tt.Output) {
				t.Errorf("ReadAll() error:unexpected feature values\ngot \t%v\nwant \t%v", res, tt.Output)
			}
		})
	}
}
