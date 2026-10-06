# Official TASSEL/rTASSEL maize tutorial data

The files in this directory originate from the Maize Genetics rTASSEL project. They are not simulated by LDblockR.

- Repository: https://github.com/maize-genetics/rTASSEL
- Pinned commit: `e8cbfcce47422319b5bfe8376cf5550ebd3b365e`
- Original directory: `inst/extdata/`
- Tutorial: https://rtassel.maizegenetics.net/articles/rTASSEL.html
- Retrieved and verified: 2026-09-02

`mdp_genotype.hmp.txt` is included byte-for-byte from the supplied file and matches the pinned upstream SHA-256 checksum. Its 281 taxa, 3,093 markers, and AGPv1 coordinates are retained. The accompanying phenotype and structure files are included for examples and are not required to read the HapMap genotype file. `PROVENANCE.json` records both the upstream and packaged byte counts and SHA-256 checksums. Two companion text files have LF line endings in the package; their tabular records are unchanged.

The HapMap allele metadata conflict with some observed calls at 48 markers. The source file is preserved without edits; filtering and genotype decoding occur during analysis.

The upstream Apache License 2.0 notice is retained in `LICENSE.rTASSEL.txt`. LDblockR does not claim ownership of these data.

Cite the original software resources:
- Bradbury et al. (2007), TASSEL. https://doi.org/10.1093/bioinformatics/btm308
- Monier et al. (2022), rTASSEL. https://doi.org/10.21105/joss.04530
- Flint-Garcia et al. (2005), maize association population. https://doi.org/10.1111/j.1365-313X.2005.02591.x
