[![PyPI version](https://badge.fury.io/py/alv.svg)](https://badge.fury.io/py/alv) 
Devel: [![Build Status](https://travis-ci.org/arvestad/alv.svg?branch=devel)](https://travis-ci.org/arvestad/alv)
# alv: a command-line alignment viewer

View you DNA or protein multiple-sequence alignments right at your command line. No need to launch a
GUI!

Note: `alv` requires Python v3.x.

Features:

* Command-line based, no GUI, so easy to script viewing of many (typically small) MSAs.
* Reads alignments in FASTA, Clustal, PHYLIP, and Stockholm formats. 
* Output is formatted to suit your terminal. You can also set the alignment width with option `-w`.
* Can color alignments of coding DNA by codon's translations to amino acids.
* Guesses sequence type (DNA/RNA/AA/coding) by default. You can override with option `-t`.
* Order sequence explicitly, alphabetically, or by sequence similarity.
* Restrict coloring to where you don't have indels or where there is a lot of conservation.

## Examples

Quick viewing of a small alignment:
```
alf msa.fa
```
This autodetects sequence type (AA, DNA, RNA, coding DNA), colors the sequences, and formats the
alignment for easy viewing in your terminal.

View three sequences, accessions `a`, `b`, and `c`, from an alignment:
```
alf -so a,b,c msa.fa
```

Feed alignment to `less`, for paging support.
```
alv -k msa.fa | less -R
```
The `-k` option ensures that `alv` keeps coloring the alignment (by default, piping
and redirection removes colors), and the `-R` option instructs `less` to interpret color codes.

## Install

Recommended installation is with `pip install alv`.

## For developers

Run `python setup.py develop test` for development install and to execute tests.

## Screenshot

### Full PFAM domain

All of the sequences in PFAM's seed alignment for PF00005

![PF00005 seed MSA](https://github.com/arvestad/alv/blob/master/doc/screenshot_PF00005.png)

### Ten peptide sequences from PF00005

![MSA from PF00005](https://github.com/arvestad/alv/blob/master/doc/screenshot_1.png)

### Seven coding DNA sequences

`alv` is autodetecting that the given DNA sequences are coding and therefore colors codons instead
of nucleotides.
![Sample screenshot](https://github.com/arvestad/alv/blob/master/doc/screenshot_2.png)
