# Default Mass Dictionary for Glycan Residues

A named numeric vector of the mass of each monosaccharide residue. Used
by default as the `mass_dict` argument in
[`calculate_mz()`](https://glycoverse.github.io/glyanno/reference/calculate_mz.md).

The names includes:

- Monosaccharides: Hex, HexNAc, dHex, dHexNAc, ddHex, Pen, HexA, HexN,
  NeuAc, NeuGc, Kdn, Neu

- Ions: H+, H, H2O, K+, Na+, NH4+, H-, Cl-, HCO3-

- Substituents: S, P

- Reducing end: red_end, which is the additional mass of the reducing
  end of the glycan.

## Usage

``` r
glyanno_mass_dict(deriv = "none", mass_type = "mono")
```

## Arguments

- deriv:

  A character scalar of the derivatization method to use. Can be "none",
  "permethyl", or "peracetyl". Default is "none".

- mass_type:

  "mono" for monoisotopic mass, "average" for average mass. Default is
  "mono".

## Value

A named numeric vector.

## Examples

``` r
glyanno_mass_dict(deriv = "none", mass_type = "mono")
#>         H       H2O        H+        K+       Na+      NH4+        H-       Cl- 
#>   1.00728  18.01056   1.00727  38.96371  22.98977  18.03382  -1.00728  34.96885 
#>     HCO3-       Hex    HexNAc      dHex   dHexNAc     ddHex       Pen      HexA 
#>  60.98014 162.05280 203.07940 146.05790 187.07980 130.05840 132.04230 176.03209 
#>      HexN     NeuAc     NeuGc       Kdn       Neu         S         P   red_end 
#> 143.05360 291.09540 307.09030 250.06890 249.08020  79.95680  79.96630  18.01056 
glyanno_mass_dict(deriv = "permethyl", mass_type = "mono")
#>         H       H2O        H+        K+       Na+      NH4+        H-       Cl- 
#>   1.00728  18.01056   1.00727  38.96371  22.98977  18.03382  -1.00728  34.96885 
#>     HCO3-       Hex    HexNAc      dHex   dHexNAc     ddHex       Pen      HexA 
#>  60.98014 204.09980 245.12630 174.08920 215.11110 144.07400 160.07360 218.07900 
#>      HexN     NeuAc     NeuGc       Kdn       Neu         S         P   red_end 
#> 217.12680 361.17370 391.18420 320.14720 333.17410  65.94120  93.98200  46.04190 
glyanno_mass_dict(deriv = "none", mass_type = "average")
#>         H       H2O        H+        K+       Na+      NH4+        H-       Cl- 
#>   1.00794  18.01524   1.00739  39.09830  22.99898  18.03850  -1.00794  35.45300 
#>     HCO3-       Hex    HexNAc      dHex   dHexNAc     ddHex       Pen      HexA 
#>  61.01680 162.14240 203.19500 146.14300 187.19480 130.14280 132.11610 176.12590 
#>      HexN     NeuAc     NeuGc       Kdn       Neu         S         P   red_end 
#> 143.14520 291.25790 307.25730 250.20530 249.21880  80.06420  79.97990  18.01524 
```
