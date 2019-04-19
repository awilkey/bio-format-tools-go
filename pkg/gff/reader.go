// Package gff reads and writes gff3 files.
// This package supports the format described in:
// http://www.sequenceontology.org/gff3.shtml
// As per the spec, gff files contain zero or more features
// of nine tab-separated fields, with the ninth column comprised of
// one or more semicolon separated fields.
//
// Feature lines that start with a # are considered comments and ignored,
// and pragma handling hasn't been implemented at this time
package gff

import (
	"bio-format-tools-go"
	"bufio"
	"bytes"
	"errors"
	"io"
	"math"
	"strconv"
	"unicode/utf8"

	//	"unicode/utf8"
	"fmt"
)

type Reader struct {
	buf        *bufio.Reader
	LineNumber uint64
	r          io.Reader
}

// All columns in a gff3 allow a "." to indicate a missing value
// This isn't and issue with columns 1-3,
const MissingValueField = math.MaxFloat64
const MissingPosField = 0
const MissingPhaseField = 3

func NewReader(r io.Reader) *Reader {
	buf := bufio.NewReader(r)
	var LineNumber uint64
	return &Reader{buf, LineNumber, r}
}

func (gr *Reader) Read() (*bio_format_tools_go.Feature, error) {
	return gr.parseFeature()
}

func (gr *Reader) parseFeature() (*bio_format_tools_go.Feature, error) {
	var line []byte
	var readErr error
	// Read next line(s), skipping comments
	for readErr == nil {
		gr.LineNumber++
		line, readErr = gr.buf.ReadBytes('\n')
		if firstRune, _ := utf8.DecodeRune(line); firstRune == '#' {
			line = nil
			continue //skip comments/pragma for now
		}
		break
	}

	// Return if read error
	if readErr != nil {
		fmt.Println(readErr, len(line))
		if len(line) == 0 && readErr == io.EOF {
			return nil, io.EOF //EOF is expected, don't bother with error
		} else if len(line) > 0 && readErr != io.EOF {
			return nil, readErr //return error
		}
	}

	fields := bytes.Split(line, []byte{'\t'})

	// Throw error if wrong number of fields
	if len(fields) != 9 {
		return nil, errors.New("wrong number of fields")
	}

	// process feature
	var feat = new(bio_format_tools_go.Feature)
	feat.Seqid = string(fields[0])
	feat.Source = string(fields[1])
	feat.Type = string(fields[2])

	if fld := string(fields[3]); fld != "." {
		feat.Start, _ = strconv.ParseUint(fld, 10, 64)
	} else {
		feat.Start = MissingPosField
	}

	if fld := string(fields[4]); fld != "." {
		feat.End, _ = strconv.ParseUint(fld, 10, 64)
	} else {
		feat.End = MissingPosField
	}

	if fld := string(fields[5]); fld != "." {
		feat.Score, _ = strconv.ParseFloat(fld, 64)
	} else {
		feat.Score = MissingValueField
	}

	if fld := string(fields[6]); fld == "+" || fld == "-" || fld == "?" || fld == "." {
		feat.Strand = fld
	} else {
		feat.Strand = "."
	}

	if fld := int8(fields[7][0]); fld > -0 && fld < 3 {
		feat.Phase = fld
	} else {
		feat.Phase = MissingPhaseField
	}

	if string(fields[8]) != "." {
		attributes := make(map[string]string)
		attrFields := bytes.Split(fields[8], []byte{';'})
		for _, attr := range attrFields {
			att := bytes.Split(attr, []byte{'='})
			if len(att) == 2 {
				attributes[string(att[0])] = string(att[1])
			}
		}
		feat.Attributes = attributes
	}

	return feat, readErr
}