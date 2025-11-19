# glyanno

The goal of glyanno is to enhance the level of information and
resolution for glycans through a hierarchical annotation approach.
Starting from a molecule mass, glyanno calculates possible glycan
compositions. Given a glycan composition, glyanno further deduces
possible glycan structures. For glycan structures with generic
monosaccharides (e.g., “Hex”, “HexNAc”), glyanno refines them into
specific types (e.g., “Glc”, “Gal”). For structures lacking linkage
information (e.g., “Gal(??-?)GalNAc(??-”), glyanno infers the most
likely linkages (e.g., “Gal(b1-3)GalNAc(a1-”).

## Installation

You can install the development version from GitHub:

``` r
remotes::install_github("glycoverse/glyanno")
```

## Documentation

- 📚 Reference:
  [Here](https://glycoverse.github.io/glyanno/reference/index.html)

## Role in `glycoverse`

Glyanno enhances the precision of glycan structures, enabling the use of
rough mass spectrometry results for detailed structural analysis. For
instance, the `glyenzy` package requires fully determined glycan
structures. Users can leverage glyanno to make educated guesses from
incomplete data, thereby enabling downstream analysis with glyenzy and
other tools that demand high-resolution structural information.

## Example

Coming soon!
