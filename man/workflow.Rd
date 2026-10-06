\name{ldblockr_workflow}
\alias{ldblockr}
\alias{ld_batch}
\alias{export_ld}
\alias{ldblockr_cli}
\title{Complete, Batch, Export, and Command-Line Workflows}
\description{Run, batch, export, or invoke the complete regional LD workflow from R or a command line.}
\usage{
ldblockr(input, output_prefix, region = NULL,
  format = c("auto", "vcf", "hapmap", "plink", "matrix"),
  samples = NULL, min_maf = 0.05, max_missing = 0.25,
  max_het = 1, min_hwe_p = 0,
  measure = c("both", "r2", "dprime"),
  r2_method = c("auto", "dosage", "haplotype"), max_distance = Inf,
  min_n = 5L,
  block_method = c("gabriel", "solid_spine", "strong", "four_gamete", "fixed", "none"),
  block_metric = c("dprime", "r2"), tag_threshold = 0.8,
  gwas = NULL, gff = NULL, special = NULL,
  output_formats = c("pdf", "svg"), plot_metric = measure,
  width = 11, height = 8.5, dpi = 300, compatible = TRUE,
  read_args = list(), ld_args = list(), block_args = list(),
  plot_args = list())

ld_batch(input, regions, output_dir, continue_on_error = TRUE, ...)

export_ld(ld, prefix, blocks = NULL, tags = NULL, compatible = TRUE)

ldblockr_cli(args = commandArgs(trailingOnly = TRUE))
}
\arguments{
\item{input}{Genotype input.}
\item{output_prefix, prefix}{Output prefix.}
\item{region, regions}{One region or multiple region strings/table.}
\item{format}{Input format.}
\item{samples}{Sample subset.}
\item{min_maf, max_missing, max_het, min_hwe_p}{Variant filters.}
\item{measure, r2_method, max_distance, min_n}{LD calculation settings; automatic r-squared uses phased haplotypes when complete.}
\item{block_method, block_metric}{Block settings.}
\item{tag_threshold}{Tag-SNP LD threshold.}
\item{gwas, gff, special}{Optional aligned tracks and labels.}
\item{output_formats, plot_metric}{Plot outputs and metrics.}
\item{width, height, dpi}{Figure settings.}
\item{compatible}{Create LDBlockShow-style tables.}
\item{read_args, ld_args, block_args, plot_args}{Additional named argument lists.}
\item{output_dir}{Batch output directory.}
\item{continue_on_error}{Continue a batch after errors.}
\item{ld}{An \code{ld_result} object.}
\item{blocks, tags}{Optional result tables.}
\item{args}{Command-line arguments.}
\item{...}{Arguments passed to \code{ldblockr()}.}
}
\value{A complete result list, batch list, named output paths, or invisible command-line result.}
\examples{
vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
\dontrun{
result <- ldblockr(vcf, "LD_demo", "chr1:1000000-1100000",
                   min_maf = 0.01, block_method = "strong",
                   output_formats = "pdf")
}
}
