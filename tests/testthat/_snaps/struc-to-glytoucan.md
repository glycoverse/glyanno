# struc_to_glytoucan rejects generic or mixed elements

    Code
      struc_to_glytoucan(generic)
    Condition
      Error in `.assert_concrete()`:
      ! `strucs` must have concrete monosaccharides (e.g. Gal, GalNAc).

---

    Code
      struc_to_glytoucan(mixed)
    Condition
      Error in `.assert_concrete()`:
      ! `strucs` must have concrete monosaccharides (e.g. Gal, GalNAc).

# struc_to_glytoucan returns NA for floating inputs

    Code
      result <- struc_to_glytoucan(floating)
    Condition
      Warning:
      `strucs` contains 1 structure with unresolved floating parts or substituents.
      i Those input structures were excluded from matching and are returned as missing values in aligned outputs.

