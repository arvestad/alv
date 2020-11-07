# Changes since last version


## v1.6.0

* New options: --only-variable and --only-variable-excluding-indels, contributed by nikostr. Constrains coloring
  to columns with variation and variation not counting indels.
* Fixed the --dotted option, which only worked with the first block for DNA sequences. Also improved the coloring
  which was too ugly in dotted mode (due to laziness).

## v1.5.0

* New option: `-d` or `--dotted`: the first sequence in the output alignment is used as a template and for positions
  in subsequent sequences that are identical, a period ('.') is output instead of a symbol.
* Adjustment: replacing blue with cyan in the DNA coloring scheme.

## v1.4.0

* New option: `-r k` or `--random-accessions k` for only showing a sample of _k_ sequences.
* New option: `-g` or `--glimpse`: display an informative cut-out of the input MSA, if it does
  not fit without scrolling or line-breaking.

## v1.3.4

* For some reason, setup.py was/is not putting proper Markdown to PyPi. Solved it halfway. Weird issue.

## v1.3.3

* Removed bug in alignmnent-type guessing: An MSA without any full codon would be guessed to be
  coding sequences.

## v1.3.2

* Ensured that alignment text color is black when starting. Helps for those other color preferences
  in the terminal.

## v1.3.1

* Adding self-promoting --cite and --methods options, noting that there is now a paper
  about `alv` in J Open Source Software.

## v1.3

* Tick marks come out as ^ on non-UTF systems, instead of an up-arrow on UTF systems (decided by looking at $LANG).
* Support for Nexus files added (thanks @SimonGreenhill)
* Added CONTRIBUTING.md, a text on how to help develop `alv`.
* Better README

## v1.2.0

* The option `-l` colors the alignment but does not break it into blocks. Suitable for piping to `less -RS`,
  as suggested by Mark McMullan <Mark.McMullan@earlham.ac.uk>.
* More indices indicated below alignments, and with an up-arrow as a tick mark.
* Added option `-sm` to allow restricting output to sequences with accessions containing a given string.

## v1.1.0

* Implemented format guessing. By default, `alv` tries to identify the alignment format. You can still override if there are parsing problems.
* You can specify what part of accessions to view with -as. This way you can avoid common prefixes and easily shorten long accessions.
* View only the first and last N characters of accessions with option -aa.
* Read from stdin with the magic file name "-".


## Introduced in version 1.0.5

* Explicitly start with "/usr/bin/env python3" instead of "/usr/bin/env python". Hoping that will
  help users that have not yet migrated to v3.
* Added v3 requirement in setup.py as well, using advice from Tobias Jakobi. Solves issue #1 on github.com.
