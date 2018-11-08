# Changes since last version

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


