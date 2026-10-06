\name{ld_blocks}
\alias{detect_ld_blocks}
\alias{select_tag_snps}
\title{Detect Haplotype Blocks and Select Tag SNPs}
\description{Call non-overlapping haplotype blocks with several definitions and select variants that greedily cover them.}
\usage{
detect_ld_blocks(ld,
  method = c("gabriel", "solid_spine", "strong", "four_gamete", "fixed", "none"),
  metric = c("dprime", "r2"), strong_cut = 0.85,
  strong_fraction = 0.9, spine_cut = 0.8, max_span = Inf,
  min_snps = 2L, fixed = NULL, gabriel_lowci = 0.7,
  gabriel_highci = 0.98, gabriel_recomb_highci = 0.9,
  gabriel_inform_fraction = 0.95, ci_step = 0.01,
  four_gamete_min_frequency = 0.01)

select_tag_snps(ld, threshold = 0.8, metric = c("r2", "dprime"),
  blocks = NULL, include_unblocked = TRUE)
}
\arguments{
\item{ld}{An \code{ld_result} object.}
\item{method}{Block definition.}
\item{metric}{Metric for applicable methods or tag selection.}
\item{strong_cut, strong_fraction}{Custom strong-pair thresholds.}
\item{spine_cut}{Solid-spine endpoint threshold.}
\item{max_span}{Maximum block length in base pairs.}
\item{min_snps}{Minimum SNPs per block.}
\item{fixed}{Fixed chromosome/start/end block table.}
\item{gabriel_lowci, gabriel_highci, gabriel_recomb_highci, gabriel_inform_fraction}{Gabriel pair and interval thresholds.}
\item{ci_step}{Likelihood grid step.}
\item{four_gamete_min_frequency}{Frequency at which a haplotype counts as observed.}
\item{threshold}{Minimum tag-SNP LD coverage.}
\item{blocks}{Optional block table for within-block selection. Both
\code{ld_blocks} output with \code{start_index}/\code{end_index} and
coordinate tables with \code{chr}/\code{start}/\code{end} are accepted.}
\item{include_unblocked}{Select unblocked variants as singleton tags.}
}
\value{A block or tag-SNP data frame.}
\examples{
vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
x <- read_vcf_region(vcf, min_maf = 0.01, quiet = TRUE)
ld <- ld_compute(x, "both")
blocks <- detect_ld_blocks(ld, method = "strong", metric = "r2",
                           strong_cut = 0.7, strong_fraction = 0.7)
select_tag_snps(ld, blocks = blocks)
}
