# deprecated db remains usable and rejects mixed filtering

    Code
      suppressWarnings(comp_to_struc("Gal(1)GalNAc(1)", db = db, glycan_type = "O-GalNAc"))
    Condition
      Error in `.resolve_annotation_db()`:
      ! Cannot combine deprecated `db` with glydb filtering arguments.
