[![PyPI version](https://badge.fury.io/py/alv.svg)](https://badge.fury.io/py/alv)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00955/status.svg)](https://doi.org/10.21105/joss.00955)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1477804.svg)](https://doi.org/10.5281/zenodo.1477804)
[![Downloads](http://pepy.tech/badge/alv)](http://pepy.tech/project/alv)
[![Build Status](https://travis-ci.org/arvestad/alv.svg?branch=master)](https://travis-ci.org/arvestad/alv)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=arvestad_alv&metric=alert_status)](https://sonarcloud.io/dashboard?id=arvestad_alv)


# alv: a command-line alignment viewer

View your DNA or protein multiple-sequence alignments right at your command line. No need to launch a
GUI!

Note: `alv` requires Python v3.4 or later. Earlier versions may also work, but this has not been
tested.

## Latest feature addition

* The command `alv -g huge_msa.fa` displays cut-out of the MSA, guaranteed to fit
  one terminal page without scrolling or MSA line breaking, that is supposed to
  give you an idea of alignment quality and contents.
* Write `alv -r 20 huge_msa.fa` to get a view of the MSA containing only 20 randomly
  selected sequences.

## Features

* Command-line based, no GUI, so easy to script viewing of many (typically small) MSAs.
* Reads alignments in FASTA, Clustal, PHYLIP, NEXUS, and Stockholm formats, from file or `stdin`.
* Output is formatted to suit your terminal. You can also set the alignment width with option `-w`.
* Can color alignments of coding DNA by codon's translations to amino acids.
* Guesses sequence type (DNA/RNA/AA/coding) by default. You can override with option `-t`.
* Order sequence explicitly, alphabetically, or by sequence similarity.
* Restrict coloring to where you don't have indels or where there is a lot of conservation.

## Install

Recommended installation is:
```
pip install --upgrade pip
pip install alv
```

If you have a half-modern BioPython installed, Python v3.4 _should_ work.
BioPython is a dependency and will only get installed automatially with `pip install alv`
if you are using Python v3.6 or later, because BioPython was apparently not on PyPi before that.


## Examples

Quick viewing of a small alignment:
```
alv msa.fa
```
This autodetects sequence type (AA, DNA, RNA, coding DNA), colors the sequences, and formats the
alignment for easy viewing in your terminal.
When applying `alv` to an alignment of coding DNA, the coding property is autodetected and colors are therefore applied to codons instead
of nucleotides.
![Seven coding DNA sequences](https://github.com/arvestad/alv/raw/master/doc/screenshot_2.png)



View three sequences, accessions `a`, `b`, and `c`, from an alignment:
```
alv -so a,b,c msa.fa
```

Feed alignment to `less`, for paging support.
```
alv -k msa.fa | less -R
```
The `-k` option ensures that `alv` keeps coloring the alignment (by default, piping
and redirection removes colors), and the `-R` option instructs `less` to interpret color codes.

## For developers

Run `python setup.py develop test` for development install and to execute tests.

## Screenshots

### Full PFAM domain

All of the sequences in PFAM's seed alignment for PF00005

![PF00005 seed MSA](https://github.com/arvestad/alv/raw/master/doc/screenshot_PF00005.png)

### Yeast sequences from PF00005

Using the option `-sm YEAST`, we reduce the alignment to the ones with a matching accession.

![Small MSA from PF00005](https://github.com/arvestad/alv/raw/master/doc/PF00005_yeast.png)
