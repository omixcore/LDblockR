# LDblockR

`LDblockR` is an installable R package for regional linkage disequilibrium (LD) analysis and visualization. It covers the main workflow of LDBlockShow while adding HapMap input, subpopulation comparisons, LD decay, haplotype frequency estimation, batch region processing, a pure R API, and automatic large-matrix rasterization.

The current submitted version is **0.0.1**. The canonical project repository, source code, and issue tracking are hosted on Gitee:
<https://gitee.com/omixcore/LDblockR>.

A public GitHub mirror is also available for code browsing and collaboration:
<https://github.com/omixcore/LDblockR>.

The package depends only on R's recommended components and does not require `ggplot2`, Bioconductor, or Java. Plain VCF/VCF.GZ files can be read directly; when `bcftools` is installed, indexed VCF region extraction is automatically accelerated.

## Main Features

| Feature | LDblockR |
|---|---|
| VCF / VCF.GZ | Supported, built-in streaming reader; automatic bcftools invocation |
| HapMap | Supported, suitable for maize population data |
| PLINK BED/BIM/FAM | Native reading |
| PLINK PED/MAP | Native reading |
| Numerical genotype matrix | Supports 0/1/2/NA |
| Subpopulation analysis | Supports arbitrary sample groups and between-group LD differences |
| LD metrics | Automatic selection of dosage or phased-haplotype \(r^2\), EM-estimated \(D'\) |
| Block detection | Gabriel, solid spine, strong-pair, four-gamete, fixed interval, no block |
| tagSNP selection | Greedy minimum-cover selection |
| Joint tracks | GWAS, GFF3/GTF, block, tagSNP, specified SNP; optional MAF track |
| Large plots | Automatic switching between vector/raster heatmaps; statistical results never subsampled |
| Graphics | Manhattan, Q-Q, regional LD; PDF, SVG, PNG, TIFF |
| Tables | Sites, pairwise LD, blocks, tags, and LDBlockShow-style files |
| Extended analysis | LD decay, haplotype frequencies, neighboring SNPs, batch regions |
| Usage modes | R functions, one-click pipeline, command line |

## Citation

If you use LDblockR in your research, please cite:

Song X, Chen P, Liang F (2026).  
*LDblockR: a cross-platform R toolkit for multi-format linkage disequilibrium analysis and integrated regional genomics visualization*.  
Frontiers in Plant Science.  
doi:10.3389/fpls.2026.1991155

BibTeX:

```bibtex
@article{Song2026LDblockR,
  author  = {Song, Xiaofeng and Chen, Pengyu and Liang, Fei},
  title   = {LDblockR: a cross-platform R toolkit for multi-format linkage disequilibrium analysis and integrated regional genomics visualization},
  journal = {Frontiers in Plant Science},
  year    = {2026},
  doi     = {10.3389/fpls.2026.1991155}
}
```

If your research is inspired by the LDBlockShow methodology, please also cite its original paper:
Dong et al. (2021), *Briefings in Bioinformatics*, DOI: 10.1093/bib/bbaa227.

## Installation

```r
# Install from GitHub
pak::pak("omixcore/LDblockR")

# Or install from Gitee
pak::pak("git::https://gitee.com/omixcore/LDblockR.git")

library(LDblockR)
packageVersion("LDblockR")  # 0.0.1
```

