# glyanno 0.4.1

## Minor improvements and bug fixes

* `com_to_struc()` now explicitly returns when the output is empty.

# glyanno 0.4.0

## Breaking changes

* `return_best = TRUE` now returns a vector with the same length as input (with `NA` for unmatched glycans) instead of a tibble with potentially fewer rows. Affected functions: `enhance_comp()`, `comp_to_struc()`, `mz_to_comp()`, `enhance_struc()`.

## Minor improvements and bug fixes

* `comp_to_struc()` now errors when concrete compositions are matched against a generic structure database, instead of silently returning empty results.
* `db` paramter cannot be of zero length now, in which case an error will be raised.

# glyanno 0.3.0

## Breaking changes

* Remove `to_level` parameter from `enhance_struc()`.

# glyanno 0.2.0

## New features

* Add `return_best` parameter to keep the most cited match when multiple compositions or structures are possible.

# glyanno 0.1.2

## Minor improvements and bug fixes

* Update dependency strategy to use the r-universe repo.

# glyanno 0.1.1

## Minor improvements and bug fixes

* Adapt to glyrepr 0.10.0.

# glyanno 0.1.0

* First GitHub release.
