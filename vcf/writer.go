package vcf

import (
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Writer allows writing gff3 files
type Writer struct {
	io.Writer
	Header bool
}

// NewWriter returns a writer after appending gff header
func NewWriter(w io.Writer) (*Writer, error) {
	return &Writer{w, false}, nil
}

func (w *Writer) WriteHeader(h Header) {

	_, _ = fmt.Fprintf(w, "##fileformat=%s\n", h.FileFormat)
	for _, val := range h.SingleVals { // Print all ##key=value lines
		_, _ = fmt.Fprintf(w, "%s\n", val)
	}
	for _, val := range h.PrintOrder { // Print all ##key=<key=value> lines
		_, _ = fmt.Fprintf(w, "%s\n", val)
	}

	_, _ = fmt.Fprintf(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(h.Genotypes) > 0 {
		gt := make([]string, len(h.Genotypes))
		for key, val := range h.Genotypes {
			gt[val] = key
		}
		_, _ = fmt.Fprintf(w, "\tFORMAT\t%s", strings.Join(gt, "\t"))
	}
}

// WriteFeature writes a single gff feature line
func (w *Writer) WriteFeature(f *Feature, h ...*Header) {
	// Write header if provided and it hasn't been printed already
	if h != nil {
		if w.Header == false {
			w.WriteHeader(*h[0])
		}
	}

	//Prep QUAL and INFO fields for pretty printing
	var qual string
	if f.Qual == MissingQualField {
		qual = "."
	} else {
		qual = strconv.FormatFloat(f.Qual, f.QualFormat, -1, 64)
	}
	info := make([]string, len(f.Info))
	for key, i := range f.InfoOrder {
		val := f.Info[key]
		if key != val {
			info[i] = fmt.Sprintf("%s=%s", key, val)
		} else {
			info[i] = fmt.Sprintf("%s", key)
		}
	}
	// print required lines
	_, _ = fmt.Fprintf(w, "\n%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s", f.Chrom, f.Pos, f.Id, f.Ref, strings.Join(f.Alt, ","), qual, f.Filter, strings.Join(info, ";"))

	// print genotype values
	if len(f.Genotypes) > 0 {
		form := make([]string, len(f.Format))
		for key, val := range f.Format {
			form[val] = key
		}
		_, _ = fmt.Fprintf(w, "\t%s\t%s", strings.Join(form, ":"), bytes.Join(f.Genotypes, []byte{'\t'}))
	}
}

// WriteAll writes all features in a slice
func (w *Writer) WriteAll(f []*Feature, h ...*Header) {
	if h[0] != nil {
		w.WriteHeader(*h[0])
	}

	for _, line := range f {
		w.WriteFeature(line, nil)
	}
}
