---
title: 'alv: a console-based viewer for molecular sequence alignments'
tags:
  - bioinformatics
  - console
  - multiple sequence alignment
  - DNA
  - RNA
  - protein
authors:
 - name: Lars Arvestad
   orcid: 0000-0001-5341-1733
   affiliation: "1, 2, 3"
affiliations:
 - name: Department of Mathematics, Stockholm University, Sweden
   index: 1
 - name: Science for Life Laboratory, Solna, Sweden
   index: 2
 - name: Swedish e-science Research Centre
   index: 3
date: 13 Aug 2018
bibliography: paper.bib
---

# Summary

The multiple sequence alignment (MSA) is a common entity in comparative analysis of molecular
sequences representing molecules such as DNA, RNA, and proteins. An MSA lines up the sequence
building blocks (letters representing nucleotides for DNA/RNA and amino acids for proteins) to form
the basis for a hypothesis of how the molecules have evolved, and is computed using, for example,
software like Clustal Omega [@sievers2014clustal], MAFFT [@katoh2013mafft], MUSCLE
[@edgar2004muscle], MACSE [@ranwez2011macse], and hmmalign [@hmmer]. MSAs have many applications,
from advanced analyses such as inferring evolutionary trees (phylogenies) or identifying function in
subsequences, to basic use like visual inspection of data. We have written a tool named `alv` to
support quick and basic viewing of MSAs [@arvestad18].

There are a number of MSA viewers available; JalView [@waterhouse2009jalview], SeaView
[@gouy2009seaview], AliView [@larsson2014aliview], and MEGA [@mega] are popular applications with many
features, including built-in analysis tools. However, due to their graphical user-interfaces, these
programs do not always work well in a command-line based workflow.  Web-based
MSA viewers are also used, for example NCBI's MSA Viewer [@ncbimsaviewer], EBI's MView [@mview], and Wasabi
[@veidenberg2015wasabi]. While offering the
advantage of not needing local software installation, yet providing analysis features,  online
tools are inconvenient when working on the command line.

Much simpler tools suffice for quick browsing of MSAs. In fact, alignment formats like PHYLIP and
Stockholm are designed to be easily read by both computers and humans, and are easily inspected with
common command-line tools (e.g., `less`) or text editors.  However, as pure text formats they lack
color, which many feel improve visual interpretation of an alignment, and suffer from a fixed
layout, which translates to suboptimal use of screen estate.

The `alv` software is an MSA viewer designed to work well in a command-line based environment and
the typical invocation is simply `alv msa.fa`. Intended use cases for `alv` includes immediate
inspection of a new alignment and quick, scriptable, browsing of many alignments. The viewer is
invoked with a straightforward command and has a number of options available. Several MSA formats
are recognized automatically (FASTA, Clustal, PHYLIP, Stockholm) and the input sequence type (DNA,
RNA, AA, or coding DNA) is guessed by default, but can also be decided when invoking `alv`.  The
output is written to stdout, with a layout adapted to the size of the current terminal and colored
to highlight similarity. For coding DNA, codons are colored according to their amino acid
translation (and several genetic codes are supported). Stop codons and frameshifts are easily
identified thanks to a highlighting color scheme.  Additional options are available to adapt the MSA
output to the user's needs.

We recommend installing `alv` using PyPi: `pip install alv`. Note that `alv` requries Python v3.2 or
later.


-![Exampel alignment of coding DNA.](screenshot.png)

# References
