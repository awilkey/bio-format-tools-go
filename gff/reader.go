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
	"bufio"
	"bytes"
	"errors"
	"io"
	"strconv"
	"unicode/utf8"
)

type Reader struct {
	buf        *bufio.Reader
	LineNumber uint64
	r          io.Reader
}

// NewReader returns a Reader.
func NewReader(r io.Reader) *Reader {
	buf := bufio.NewReader(r)
	var LineNumber uint64
	return &Reader{buf, LineNumber, r}
}

// Read returns a pointer to a Feature. Input is assumed to be a properly formed gff3
func (gr *Reader) Read() (*Feature, error) {
	return gr.parseFeature()
}

// ReadAll returns a slice of pointers to Features from an input of one-or-more lines
func (gr *Reader) ReadAll() (features []*Feature, err error) {
	for {
		feature, err := gr.parseFeature()
		if feature != nil {
			feat := feature
			features = append(features, feat)
		}
		if err == io.EOF {
			return features, err
		}
		if err != nil {
			return features, err
		}
	}
}

//func (gr *Reader) validateFeature(fields [][]byte) error{
// TODO: Add ReadAndValidate
//	if field := strings.TrimSpace(string(fields[0])); field == "" || field == "." {
//		return errors.New("column 1 (seqid) must be defined")
//	}
//
//	if field := strings.TrimSpace(string(fields[1])); field == "" {
//		return errors.New("column 2 (source) must not be empty")
//	}
//
//	if field := strings.TrimSpace(string(fields[2])); field == "" || field == "." {
//		return errors.New("column 3 (type) must be defined")
//	}
//
//	if field := strings.TrimSpace(string(fields[2])); field == "" || field == "." {
//		return errors.New("column 3 (type) must be defined")
//	}
//
//	if field := strings.TrimSpace(string(fields[3])); field == "." || field == "" {
//		return errors.New("column 4 (start) must be defined")
//	}
//
//	start,_ := strconv.ParseUint(string(fields[3]), 10, 64)
//	if start < 1 {
//		return errors.New("column 4 (start) must be one-based")
//	}
//
//	if  field := strings.TrimSpace(string(fields[4])); field == "." || field == "" {
//		return errors.New("column 5 (end) must be defined")
//	}
//
//	end,_ := strconv.ParseUint(string(fields[4]), 10, 64)
//	if end < 1 {
//		return errors.New("column 5 (end) must be one-based")
//	}
//
//	if start > end {
//		return errors.New("column 5 (end) must be greater than or equal to column 4 (start)")
//	}
//
//	return nil
//}

func (gr *Reader) parseFeature() (*Feature, error) {
	var line []byte
	var readErr error
	// Read next line(s), skipping comments
	for readErr == nil {
		gr.LineNumber++
		line, readErr = gr.buf.ReadBytes('\n')
		if firstRune, _ := utf8.DecodeRune(line); firstRune == '#' || bytes.TrimSpace(line) == nil {
			line = nil
			continue //skip comments/pragma for now
		}

		break
	}

	// Return if read error
	if readErr != nil {
		if len(line) == 0 && readErr == io.EOF {
			return nil, io.EOF //EOF is expected, don't bother with error
		} else if len(line) > 0 && readErr != io.EOF {
			return nil, readErr //return error
		}
	}

	fields := bytes.Split(line, []byte{'\t'})

	// Throw error if wrong number of fields
	if !(len(fields) == 9 || len(fields) == 8) {
		return nil, errors.New("wrong number of fields")
	}

	// process feature
	var feat = new(Feature)
	feat.Seqid = string(fields[0])
	feat.Source = string(fields[1])
	feat.Type = string(fields[2])
	feat.Start, _ = strconv.ParseUint(string(fields[3]), 10, 64)
	feat.End, _ = strconv.ParseUint(string(fields[4]), 10, 64)
	if fld := string(fields[5]); fld != "." {
		feat.Score, _ = strconv.ParseFloat(fld, 64)
	} else {
		feat.Score = MissingScoreField
	}
	feat.Strand = string(fields[6])
	if fld, _ := strconv.ParseInt(string(fields[7]), 10, 8); string(fields[7]) != "." && fld >= 0 && fld < 3 {
		feat.Phase = int8(fld)
	} else {
		feat.Phase = MissingPhaseField
	}

	if len(fields) == 9 {
		attributes := map[string]string{}
		if string(fields[8]) != "." {
			attrFields := bytes.Split(fields[8], []byte{';'})
			for _, attr := range attrFields {
				att := bytes.Split(attr, []byte{'='})
				if len(att) == 2 {
					attributes[string(bytes.TrimSpace(att[0]))] = string(bytes.TrimSpace(att[1])) //Clean leading and trailing whitespace
				}
			}
		}
		feat.Attributes = attributes
	}

	return feat, readErr
}
