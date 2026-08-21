# enhance_struc rejects floating inputs and database candidates

    Code
      best <- enhance_struc(input, db = db, return_best = TRUE)
    Condition
      Warning:
      `strucs` contains 1 structure with unresolved floating parts or substituents.
      i Those input structures were excluded from matching and are returned as missing values in aligned outputs.
      Warning:
      `db` contains 1 structure with unresolved floating parts or substituents.
      i Those database structures were excluded from matching.

# enhance_struc_denovo rejects a non-concrete fallback database

    Code
      enhance_struc_denovo(input, fallback_db = input)
    Condition
      Error in `.prepare_denovo_struc_db()`:
      ! `fallback_db` must contain only concrete structures.

# enhance_struc_denovo requires generic topological inputs

    Code
      enhance_struc_denovo(concrete)
    Condition
      Error in `.assert_denovo_strucs()`:
      ! `strucs` must contain only generic topological structures.
      i Missing values are allowed, but every other element must have generic monosaccharides and no linkage information.

---

    Code
      enhance_struc_denovo(mixed)
    Condition
      Error in `.assert_denovo_strucs()`:
      ! `strucs` must contain only generic topological structures.
      i Missing values are allowed, but every other element must have generic monosaccharides and no linkage information.

---

    Code
      enhance_struc_denovo(partial)
    Condition
      Error in `.assert_denovo_strucs()`:
      ! `strucs` must contain only generic topological structures.
      i Missing values are allowed, but every other element must have generic monosaccharides and no linkage information.

---

    Code
      enhance_struc_denovo(intact)
    Condition
      Error in `.assert_denovo_strucs()`:
      ! `strucs` must contain only generic topological structures.
      i Missing values are allowed, but every other element must have generic monosaccharides and no linkage information.

# enhance_struc_denovo rejects floating inputs

    Code
      result <- enhance_struc_denovo(c("Hex(??-?)HexNAc(??-", floating), fallback_db = fallback_db)
    Condition
      Warning:
      `strucs` contains 1 structure with unresolved floating parts or substituents.
      i Those input structures were excluded from matching and are returned as missing values in aligned outputs.

