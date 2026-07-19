# enhance_struc de-novo requires a supported target level

    Code
      enhance_struc(input, method = "denovo")
    Condition
      Error in `enhance_struc()`:
      ! `to_level` must be provided when `method` is "denovo".

---

    Code
      enhance_struc(input, method = "denovo", to_level = "intact")
    Condition
      Error in `enhance_struc()`:
      ! De-novo enhancement to "intact" is not yet supported.

