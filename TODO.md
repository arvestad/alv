# Planned fixes and features

* Better screenshots
* Review how colors work in different terminals. 
    - Need to differentiate between dark and bright backgrounds.
    - Accession color should be explicit
    - Differentiate between coloring of sequences and "context" (accessions, indices, etc).
* Add support for restricting to a sub-alignment.
* Add option --glimpse.
* Explicitly choose parts of an alignment to view/color.
* Pypi.org does not handle MarkDown in the README. Should look for a solution. 

# Considered features

* Automatic paging using (e.g.) `less -R` _if_ alignment does not fit in one screen and env variable `ALV_PAGER` is set.
* Coloring of motifs
* Option -st, --sort-by-tree: Sort sequences based on a simply estimated tree, à la Belvu. Should be easy to implement using Dendropy.
