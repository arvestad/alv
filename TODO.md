# Planned fixes and features

* Review how colors work in different terminals. The color scheme working for me looked strange when Vilde was
  showing me her work.
* Make it possible to color alignments without breaking alignment into blocks. I.e., a "pure"
  version of `alv -k -w 99999 msa.fa | less -SR`. Suggested by Mark McMullan <Mark.McMullan@earlham.ac.uk>.
* Add support for restricting to a sub-alignment.
* Add option --glimpse.
* Explicitly choose parts of an alignment to view/color.


# Considered features

* Automatic paging using (e.g.) `less -R` _if_ alignment does not fit in one screen and env variable `ALV_PAGER` is set.
* Coloring of motifs
