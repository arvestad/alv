# alv: a command-line alignment viewer

View you DNA or protein multiple-sequence alignments right at your command line. No need to launch a
GUI!

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