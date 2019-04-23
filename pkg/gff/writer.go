package gff

import (
	"fmt"
	"io"
)

// Writer allows writing gff3 files
type Writer struct {
	io.Writer
}

// NewWriter returns a writer after appending gff header
func NewWriter(w io.Writer) (*Writer, error) {
	_, _ = fmt.Fprintf(w, "##gff-version 3.2.1\n")
	return &Writer{w}, nil
}

// WriteFeature writes a single gff feature line
func (w *Writer) WriteFeature(f *Feature) {
	_, _ = fmt.Fprintln(w, f)
}

// WriteAll writes all features in a slice
func (w *Writer) WriteAll(f []*Feature) {
	for _, line := range f {
		_, _ = fmt.Fprintln(w, line)
	}
}
