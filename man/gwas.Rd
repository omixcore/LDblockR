\name{gwas_mlm}
\alias{gwas_mlm}
\alias{print.gwas_mlm_result}
\alias{gwas_ld_region}
\alias{print.gwas_region}
\alias{read_tassel_phenotype}
\alias{read_tassel_kinship}
\alias{example_data}
\title{Fast mixed linear model GWAS and built-in demonstration data}
\description{
Fits a single-variance-component mixed linear model using a P3D-style null
fit. The kinship eigen-decomposition and variance ratio are reused for all
markers, and marker tests are evaluated in matrix chunks. The returned table
contains marker PVE and model PVE.
}
\usage{
gwas_mlm(genotype, phenotype, trait = NULL, kinship = NULL,
  covariates = NULL, n_pc = 0L, map = NULL, sample_ids = NULL,
  replicate = c("expand", "mean"), min_maf = 0.05,
  max_missing = 0.2, chunk_size = 256L, max_iter = 60L,
  verbose = TRUE)
\method{print}{gwas_mlm_result}(x, ...)
read_tassel_phenotype(file)
read_tassel_kinship(file, samples = NULL)
example_data(name = "tassel")
gwas_ld_region(gwas, cutline = NULL, flank = 250000, min_width = 0,
  max_width = Inf, chromosome = NULL)
}
\arguments{
  \item{genotype}{An \code{ld_data} object, 0/1/2 matrix, HapMap file, or
  RDS file.}
  \item{phenotype}{A data frame, TASSEL phenotype file, or named numeric
  vector.}
  \item{trait}{Trait column; the first numeric column is used when omitted.}
  \item{kinship}{Optional matrix or TASSEL kinship file.}
  \item{covariates}{Fixed-effect columns, a matrix, or \code{NULL} to use
  \code{location}, \code{Q1}, \code{Q2}, and \code{Q3} when present.}
  \item{n_pc}{Number of genotype PCs appended to the fixed effects.}
  \item{map}{Optional marker map for a matrix input.}
  \item{sample_ids}{Optional matrix-row sample IDs.}
  \item{replicate}{Keep repeated phenotype rows or average by taxon.}
  \item{min_maf}{Minimum minor allele frequency.}
  \item{max_missing}{Maximum marker missing-call fraction.}
  \item{chunk_size}{Number of markers per vectorized scan chunk.}
  \item{max_iter}{Compatibility argument retained for the variance-profile API.}
  \item{verbose}{Print null-model and progress summaries.}
  \item{x}{A result returned by \code{gwas_mlm}.}
  \item{file}{A TASSEL phenotype or kinship text file.}
  \item{samples}{Optional IDs to retain from a kinship matrix.}
  \item{name}{Built-in data bundle name: \code{"tassel"} or \code{"regional"}.}
  \item{gwas}{A GWAS table or file accepted by \code{read_gwas()}.}
  \item{cutline}{The \code{-log10(P)} threshold; \code{NULL} uses the
  Bonferroni threshold.}
  \item{flank}{Bases added on each side of the selected hits.}
  \item{min_width}{Minimum returned interval width in bases.}
  \item{max_width}{Optional maximum interval width.}
  \item{chromosome}{Optional chromosome to force.}
  \item{...}{Unused.}
}
\value{
\code{gwas_mlm} returns a data frame with class \code{gwas_mlm_result},
including \code{chr}, \code{pos}, \code{id}, \code{p}, \code{logp},
\code{beta}, \code{se}, \code{PVE}, and \code{model_PVE}. Model metadata are
stored in \code{attr(result, "model")}. The readers return normalized tables
or a named symmetric matrix; \code{example_data} returns installed file paths.
The \code{"regional"} bundle includes the VCF/HapMap, regional GWAS, GFF3,
fixed-block, and highlighted-SNP files used by the plotting examples.
\code{gwas_ld_region()} returns a \code{gwas_region} list with a region string,
lead marker, threshold-hit IDs, and the selected GWAS table. If no point exceeds
the threshold, the most significant marker is used as a fallback. Distant
threshold peaks on the same chromosome are split into local clusters so that a
single LD panel does not span unrelated loci.
}
\details{
The tested model is \eqn{y = X\beta + Zu + e}, with
\eqn{u \sim N(0, K\sigma_g^2)}. The null variance ratio
\eqn{\delta = \sigma_e^2/\sigma_g^2} is profiled once. For each marker the
whitened marker is residualized against \code{X}; the association test is a
one-degree-of-freedom GLS t test. Marker \code{PVE} is the partial reduction
in the weighted residual sum of squares, while \code{model_PVE} is
\eqn{1/(1+\delta)}. This is the same P3D principle as GAPIT, implemented with
an independent matrix-chunked scan.
}
\examples{
\dontrun{
paths <- example_data()
pheno <- read_tassel_phenotype(paths[["mdp_phenotype"]])
gwas <- gwas_mlm(paths[["mdp_genotype"]], pheno, trait = "EarHT",
  covariates = c("location", "Q1", "Q2", "Q3"),
  replicate = "expand", min_maf = 0.05, chunk_size = 256)
head(gwas[order(gwas$p), ])
attr(gwas, "model")$model_PVE
region <- gwas_ld_region(gwas, cutline = -log10(0.05 / nrow(gwas)))
region$region
}
}
