\name{ld_analysis}
\alias{ld_compute}
\alias{ld_neighbors}
\alias{ld_by_group}
\alias{ld_compare}
\alias{haplotype_frequencies}
\alias{ld_decay}
\alias{plot_ld_decay}
\title{Calculate, Query, Compare, and Summarize Linkage Disequilibrium}
\description{Calculate pairwise LD and derive neighbors, group contrasts, haplotype frequencies, and distance-decay summaries.}
\usage{
ld_compute(x, measure = c("r2", "dprime", "both"),
  r2_method = c("auto", "dosage", "haplotype"), max_distance = Inf,
  min_n = 5L, ci = FALSE, ci_step = 0.01, use_phased = TRUE)

ld_neighbors(ld, variant, threshold = 0.8, metric = c("r2", "dprime"))

ld_by_group(x, groups, ...)

ld_compare(grouped, metric = c("r2", "dprime"), reference = 1L)

haplotype_frequencies(x, variants = NULL, min_frequency = 0.01,
  drop_missing = TRUE)

ld_decay(ld, metric = c("r2", "dprime"), n_bins = 30L,
  bin_width = NULL, max_distance = NULL)

plot_ld_decay(x, col = "#B2182B", ribbon_col = "#FDD0A2",
  xlab = "Physical distance (bp)", ylab = NULL, ...)
}
\arguments{
\item{x}{An \code{ld_data} object; for plotting, an \code{ld_decay} object.}
\item{measure}{Statistics to calculate.}
\item{r2_method}{Automatic selection, dosage correlation, or phased-haplotype correlation. Automatic mode uses haplotypes when phase is complete.}
\item{max_distance}{Maximum physical pair distance.}
\item{min_n}{Minimum pairwise sample count.}
\item{ci}{Calculate profile-likelihood D-prime intervals.}
\item{ci_step}{D-prime likelihood grid step.}
\item{use_phased}{Use complete phase information when available.}
\item{ld}{An \code{ld_result} object.}
\item{variant}{Focal variant ID, index, or position.}
\item{threshold}{Minimum LD value.}
\item{metric}{LD metric.}
\item{groups}{Named vector, ordered vector, or sample/group table.}
\item{grouped}{Result from \code{ld_by_group()}.}
\item{reference}{Reference group name or index.}
\item{variants}{Variant IDs or indices.}
\item{min_frequency}{Minimum reported haplotype frequency.}
\item{drop_missing}{Drop incomplete haplotypes.}
\item{n_bins}{Number of distance bins.}
\item{bin_width}{Optional fixed bin width.}
\item{col, ribbon_col}{Decay line and ribbon colors.}
\item{xlab, ylab}{Axis labels.}
\item{...}{Additional arguments passed to \code{ld_compute()} or graphics.}
}
\value{An \code{ld_result}, query table, grouped list, comparison list, haplotype table, decay table, or invisible plot input, as appropriate.}
\examples{
vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
x <- read_vcf_region(vcf, min_maf = 0.01, quiet = TRUE)
ld <- ld_compute(x, "both")
head(ld_neighbors(ld, "snp21", threshold = 0.5))
head(ld_decay(ld))
}
