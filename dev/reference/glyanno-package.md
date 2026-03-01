# glyanno: Hierarchical Glycan Annotation

Hierarchical annotation of glycans based on molecule mass, composition,
and structure. Starting from a molecule mass, glyanno calculates
possible glycan compositions. Given a glycan composition, glyanno
further deduces possible glycan structures. For glycan structures with
generic monosaccharides (e.g., "Hex", "HexNAc"), glyanno refines them
into specific types (e.g., "Glc", "Gal"). For structures lacking linkage
information (e.g., "Gal(??-?)GalNAc(??-"), glyanno infers the most
likely linkages (e.g., "Gal(b1-3)GalNAc(a1-").

## See also

Useful links:

- <https://glycoverse.github.io/glyanno/>

## Author

**Maintainer**: Bin Fu <23110220018@m.fudan.edu.cn>
([ORCID](https://orcid.org/0000-0001-8567-2997))
