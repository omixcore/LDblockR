\name{LDblockR-package}
\alias{LDblockR}
\alias{LDblockR-package}
\title{LDblockR: Regional Linkage Disequilibrium Analysis and Visualization}
\description{
Read VCF, HapMap, PLINK, or matrix genotype data; calculate r-squared and
Lewontin's D-prime; identify haplotype blocks and tag variants; fit a fast
P3D-principle mixed linear model GWAS with PVE; and draw aligned LD,
association, gene, block, and tag-SNP tracks.
}
\details{
The package has no mandatory non-base R dependencies. See \code{gwas_mlm()} for
the cached-kinship MLM scan, \code{ldblockr()} for the one-call LD workflow,
and \code{plot_ld()} for publication-ready figures.
}
\docType{package}
\author{Xiaofeng Song, Pengyu Chen, and Fei Liang}
\references{
Dong SS et al. (2021) LDBlockShow: a fast and convenient tool for visualizing
linkage disequilibrium and haplotype blocks based on variant call format files.
\emph{Briefings in Bioinformatics}. doi:10.1093/bib/bbaa227.
}
