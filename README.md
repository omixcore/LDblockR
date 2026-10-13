# LDblockR

A pure R toolkit for regional linkage disequilibrium analysis and integrated genomics visualization. LDblockR supports VCF/VCF.GZ, HapMap, PLINK, and numeric genotype matrices on Windows, macOS, and Linux.

[**Tutorial**](docs/TUTORIAL.md) · [Validation](inst/VALIDATION.md) · [GWAS validation](inst/GWAS_VALIDATION.md) · [Changelog](NEWS.md)

## Features

- Calculate dosage or phased-haplotype r² and Lewontin's D′.
- Detect LD blocks, select tag SNPs, and summarize LD decay and haplotype frequencies.
- Compare sample groups and process multiple regions.
- Combine regional GWAS, gene models, blocks, and LD heatmaps; export PDF, SVG, PNG, or TIFF.
- Run mixed-model GWAS and generate Manhattan and Q-Q plots through an R API or command line.

## Installation

Requires R >= 4.1.0; no compiler is needed.

```r
# Install pak first if needed: install.packages("pak")
pak::pak("omixcore/LDblockR")
library(LDblockR)
```

[Local installation and Gitee alternative](docs/TUTORIAL.md#installation).

## Quick Start

```r
library(LDblockR)
regional <- example_data("regional")

result <- ldblockr(
  input = regional[["vcf"]],
  output_prefix = "LDblockR_demo",
  region = "chr1:1000000-1100000",
  min_maf = 0.01,
  measure = "both",
  block_method = "gabriel",
  gwas = regional[["regional_gwas"]],
  gff = regional[["gff3"]],
  output_formats = c("pdf", "svg", "png")
)
result$files
```

See the [full tutorial](docs/TUTORIAL.md) for input formats, plotting options, TASSEL examples, and reproducible workflows.

## Project

[GitHub](https://github.com/omixcore/LDblockR) · [Gitee](https://gitee.com/omixcore/LDblockR) · [Issues](https://github.com/omixcore/LDblockR/issues)

## Citation

If you use LDblockR in your research, please cite:

> Song, X., Chen, P., and Liang, F. (2026). LDblockR: a cross-platform R toolkit for multi-format linkage disequilibrium analysis and integrated regional genomics visualization. *Frontiers in Plant Science*. DOI: [10.3389/fpls.2026.1991155](https://doi.org/10.3389/fpls.2026.1991155).

[Article page](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1991155/full).

```bibtex
@article{Song2026LDblockR,
  author  = {Song, Xiaofeng and Chen, Pengyu and Liang, Fei},
  title   = {{LDblockR}: a cross-platform {R} toolkit for multi-format linkage disequilibrium analysis and integrated regional genomics visualization},
  journal = {Frontiers in Plant Science},
  year    = {2026},
  doi     = {10.3389/fpls.2026.1991155},
  url     = {https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1991155/full}
}
```

If your research is inspired by the LDBlockShow methodology, please also cite its original paper:
Dong et al. (2021), *Briefings in Bioinformatics*, DOI: [10.1093/bib/bbaa227](https://doi.org/10.1093/bib/bbaa227).
