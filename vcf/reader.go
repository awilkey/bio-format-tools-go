// Package gff reads and writes gff3 files.
// This package supports the format described in:
// http://www.sequenceontology.org/gff3.shtml
// As per the spec, gff files contain zero or more features
// of nine tab-separated fields, with the ninth column comprised of
// one or more semicolon separated fields.
//
// Feature lines that start with a # are considered comments and ignored,
// and pragma handling hasn't been implemented at this time
package vcf

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
)

type Reader struct {
	buf        *bufio.Reader
	Header     *Header
	LineNumber uint64
	r          io.Reader
}

// NewReader returns a Reader.
func NewReader(r io.Reader) (*Reader, error) {
	buf := bufio.NewReader(r)
	var LineNumber uint64
	var line []byte
	var readErr error
	foundHeader := false
	h := NewHeader()
	//read header/meta
	for readErr == nil {
		LineNumber++
		line, readErr = buf.ReadBytes('\n')
		if LineNumber == 1 {
			metaFields, _, _, err := parseLineToMeta(line)
			if err != nil || metaFields["FieldType"] != "fileformat" {
				readErr = errors.New("fileformat is not vcf")
				break
			} else {
				h.FileFormat = metaFields["ID"]
			}
		} else if bytes.HasPrefix(line, []byte("##")) { // Meta directive
			metaFields, metaOrder, formatted, err := parseLineToMeta(line)
			if err != nil {
				readErr = err
				break
			} else if LineNumber == 1 {
				if metaFields["FieldType"] != "fileformat" {
					readErr = errors.New("fileformat is not vcf")
					break
				} else {
					h.FileFormat = metaFields["ID"]
				}
			} else {
				if formatted == false { // Meta directive is simple ##key=value
					meta := SingleValMeta{FieldType: metaFields["FieldType"], Id: metaFields["ID"]}
					h.SingleVals = append(h.SingleVals, &meta)
				} else { // Meta directive is ##key=<key=value,...>
					var meta Meta
					meta.FieldOrder = metaOrder[1:]

					for _, field := range metaOrder { // load meta struct
						switch field {
						case "FieldType":
							meta.FieldType = metaFields[field]
						case "ID":
							meta.Id = metaFields[field]
						case "Number":
							meta.Number = metaFields[field]
						case "Type":
							meta.Type = metaFields[field]
						case "Description":
							meta.Description = metaFields[field]
						case "URL":
							meta.Url = metaFields[field]
						default:
							if meta.Optional == nil {
								meta.Optional = make(map[string]string)
							}
							meta.Optional[field] = metaFields[field]
						}
					}

					switch metaFields["FieldType"] {
					case "META":
						h.Metas = append(h.Metas, &meta)
					case "INFO":
						h.Infos = append(h.Infos, &meta)
					case "FILTER":
						h.Filters = append(h.Filters, &meta)
					case "FORMAT":
						h.Formats = append(h.Formats, &meta)
					case "ALT":
						h.Alts = append(h.Alts, &meta)
					case "SAMPLE":
						h.Samples = append(h.Samples, &meta)
					case "assembly":
						h.Assemblies = append(h.Assemblies, &meta)
					case "contig":
						h.Contigs = append(h.Contigs, &meta)
					case "pedigree":
						h.Pedigrees = append(h.Pedigrees, &meta)
					default:
						h.Others = append(h.Others, &meta)
					}

					h.PrintOrder = append(h.PrintOrder, &meta)

				}
			}
		} else if bytes.HasPrefix(line, []byte("#")) { //header
			foundHeader = true
			header := bytes.Split(line, []byte("\t"))
			if len(header) < 8 {
				readErr = errors.New("header has too few columns to be minimum vcf")
			} else if len(header) == 9 { //either no genotypes or filter + at least one genotype
				if string(header[8]) != "FORMAT" {
					readErr = errors.New("header needs a FORMAT column before adding genotypes")
				} else {
					readErr = errors.New("header FORMAT column must be followed by at least one genotype")
				}
			} else if len(header) >= 10 {
				h.Genotypes = make(map[string]uint64, len(header)-9)
				for i, genotype := range header[9:] {
					h.Genotypes[string(genotype)] = uint64(i)
				}
			}
			break
		} else { //not header or meta, can start parsing actual file
			break
		}
	}
	if !foundHeader {
		readErr = errors.New("no header line present")
	}

	// error reading header
	if readErr != nil {
		if len(line) == 0 && readErr == io.EOF {
			return nil, io.EOF //EOF is expected, don't bother with error
		} else if len(line) > 0 && readErr != io.EOF {
			return nil, readErr //return error
		}
	}

	return &Reader{buf, h, LineNumber, r}, nil
}

func parseLineToMeta(meta []byte) (map[string]string, []string, bool, error) {
	meta = bytes.TrimLeft(bytes.TrimSpace(meta), "#")
	line := bytes.SplitN(meta, []byte("="), 2)
	var formatted bool
	metaValues := make(map[string]string)
	var metaOrder []string

	if len(line) == 1 {
		return nil, nil, false, errors.New("invalid meta header")
	} else {
		if !bytes.HasPrefix(line[1], []byte("<")) {
			formatted = false
			metaValues["FieldType"] = string(line[0])
			metaValues["ID"] = string(line[1])
			metaOrder = append(metaOrder, "FieldType", "ID")
		} else {
			formatted = true
			metaValues["FieldType"] = string(bytes.TrimLeft(line[0], "#"))
			metaOrder = append(metaOrder, "FieldType")
			fields := bytes.Split(bytes.Trim(bytes.TrimSpace(line[1]), "<>"), []byte(","))
			for _, field := range fields {
				info := bytes.Split(field, []byte("="))
				if len(info) == 1 {
					continue
				} else {
					metaValues[string(info[0])] = string(info[1])
					metaOrder = append(metaOrder, string(info[0]))
				}
			}
		}
	}

	return metaValues, metaOrder, formatted, nil
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

// parseFeature from a VCF line
func (gr *Reader) parseFeature() (*Feature, error) {
	var line []byte
	var readErr error

	gr.LineNumber++
	line, readErr = gr.buf.ReadBytes('\n')

	// Return if read error
	if readErr != nil {
		if len(line) == 0 && readErr == io.EOF {
			return nil, io.EOF //EOF is expected, don't bother with error
		} else if len(line) > 0 && readErr != io.EOF {
			return nil, readErr //return error
		}
	}

	fields := bytes.Split(line, []byte{'\t'})

	if flen := len(fields); flen < 8 || flen == 9 || (flen >= 10 && flen != (len(gr.Header.Genotypes)+1+8)) { // Error if not enough fields in line
		err := fmt.Sprintf("too few columns in feature line")
		return nil, errors.New(err)
	}

	var feat Feature
	for i := range fields { // Trim extra leading/trailing spaces for each field
		fields[i] = bytes.TrimSpace(fields[i])
	}
	// Populate required fields
	feat.Chrom = string(fields[0])
	feat.Pos, _ = strconv.ParseUint(string(fields[1]), 10, 64)
	feat.Id = string(fields[2])
	feat.Ref = string(fields[3])

	alt := bytes.Split(fields[4], []byte{','})
	feat.Alt = make([]string, len(alt))
	for i := range alt {
		feat.Alt[i] = string(alt[i])
	}

	feat.Qual, _ = strconv.ParseFloat(string(fields[5]), 64)
	if bytes.IndexAny(fields[5], "eE") != -1 {
		feat.QualFormat = byte('e')
	} else {
		feat.QualFormat = byte('f')
	}

	feat.Filter = string(fields[6])

	infos := bytes.Split(fields[7], []byte{';'})
	feat.Info = make(map[string]string, len(infos))
	feat.InfoOrder = make(map[string]int, len(infos))
	for i := range infos {
		bytes.TrimSpace(infos[i])
		curInf := bytes.Split(infos[i], []byte{'='})
		if len(curInf) == 1 {
			feat.Info[string(curInf[0])] = string(curInf[0])
		} else {
			feat.Info[string(curInf[0])] = string(curInf[1])
		}
		feat.InfoOrder[string(curInf[0])] = i
	}

	if len(fields) > 8 { // if more than eight fields, populate genotype
		fmts := bytes.Split(fields[8], []byte{':'})
		feat.Format = make(map[string]int, len(fmts))
		for i := range fmts {
			feat.Format[string(fmts[i])] = i
		}

		feat.Genotypes = fields[9:]
	}

	return &feat, readErr
}
