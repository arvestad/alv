# Planned fixes and features

* Remove the `args` argument to Alignment.blocks().
* Add support for restricting to a sub-alignment.
* Add option --glimpse.
* Implement format guessing and make it the default.

# Considered features

* Automatic paging using (e.g.) `less -R` _if_ alignment does not fit in one screen and env variable `ALV_PAGER` is set.