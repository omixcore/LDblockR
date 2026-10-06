# Core LD validation record

This record documents the core LD benchmark used for the manuscript. The
final 0.0.1 release should be run through `R CMD check` before public release.

Validation date: 2026-08-24

Runtime used for release validation:

- R 4.5.1
- Ubuntu 24.04, x86_64
- LDblockShow 1.41 from the official `hewm2008/LDBlockShow` repository

## Package checks

- Clean source installation: passed.
- `R CMD check --no-manual`: 0 errors and 0 warnings.
- Base-R integration tests: passed.
- PDF, SVG, PNG, and TIFF devices: passed.
- Manhattan and Q-Q GWAS diagnostic devices: passed.
- Automatic raster rendering of a 260-variant heatmap: passed.

## Tested input paths

- Phased VCF.GZ regional reader.
- HapMap reader with IUPAC and diploid allele calls.
- PLINK SNP-major BED/BIM/FAM decoder, including missing and heterozygous calls.
- Samples-by-variants 0/1/2 matrix plus map.
- GFF3 gene/transcript/CDS/UTR annotation.
- Headered association table.

## Tested analyses

- Dosage and phased-haplotype r-squared.
- EM-estimated D-prime and profile-likelihood confidence intervals.
- Gabriel, solid-spine, custom strong-pair, four-gamete, and fixed blocks.
- Greedy tag-SNP coverage.
- Sample-group LD matrices and group-minus-reference comparisons.
- LD neighbors, decay summaries, and haplotype frequencies.
- Modern tables and native triangular LDBlockShow-compatible table exports.
- Complete one-call workflow.

## Numerical comparison with official LDBlockShow

The bundled phased example VCF contains 60 diploid samples and 42 variants.
Both programs analyzed all 861 unique variant pairs with MAF 0.01. LDblockShow
reported values to three decimal places.

| Statistic | Maximum absolute difference | Mean absolute difference | Pairs within 0.00051 |
|---|---:|---:|---:|
| phased-haplotype r-squared | 0.0004993138 | 0.0002348542 | 100% |
| D-prime | 0.0005000000 | 0.0002371987 | 100% |

Thus every LDblockR value agreed with the corresponding official LDblockShow
value to the precision written by LDblockShow.
