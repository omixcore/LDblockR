\name{genotype_io}
\alias{parse_region}
\alias{as_ld_data}
\alias{read_vcf_region}
\alias{read_hapmap}
\alias{read_plink}
\alias{read_genotypes}
\alias{filter_variants}
\alias{subset_ld_data}
\title{Read, Construct, Filter, and Subset Genotype Data}
\description{Import common population-genotype formats or construct, filter, and subset an in-memory genotype object.}
\usage{
parse_region(region)

as_ld_data(genotypes, map = NULL, sample_ids = NULL,
  haplotypes = NULL, source = NULL)

read_vcf_region(file, region = NULL, samples = NULL, min_maf = 0.05,
  max_missing = 0.25, max_het = 1, min_hwe_p = 0,
  pass_only = FALSE, biallelic_snps_only = TRUE,
  backend = c("auto", "base", "bcftools"), max_variants = Inf,
  on_excess = c("error", "thin"), quiet = FALSE)

read_hapmap(file, region = NULL, samples = NULL, min_maf = 0.05,
  max_missing = 0.25, max_het = 1, min_hwe_p = 0, quiet = FALSE)

read_plink(prefix, region = NULL, samples = NULL, min_maf = 0.05,
  max_missing = 0.25, max_het = 1, min_hwe_p = 0, quiet = FALSE)

read_genotypes(x, format = c("auto", "vcf", "hapmap", "plink", "matrix"), ...)

filter_variants(x, min_maf = 0.05, max_maf = 0.5,
  max_missing = 0.25, max_het = 1, min_hwe_p = 0, quiet = FALSE)

subset_ld_data(x, variants = NULL, samples = NULL)
}
\arguments{
\item{region}{A region string in \code{chr:start-end} or \code{chr:start:end} form.}
\item{genotypes}{Samples-by-variants numeric matrix coded 0, 1, 2, or \code{NA}.}
\item{map}{Variant table or numeric position vector.}
\item{sample_ids}{Optional sample identifiers.}
\item{haplotypes}{Optional zero/one matrix with two haplotype rows per sample.}
\item{source}{Optional source description.}
\item{file}{VCF or HapMap file, optionally gzip-compressed.}
\item{samples}{Sample IDs or a one-column sample-list file.}
\item{min_maf, max_maf}{Minor allele-frequency bounds.}
\item{max_missing}{Maximum missing-call fraction.}
\item{max_het}{Maximum heterozygous-call fraction.}
\item{min_hwe_p}{Minimum chi-square Hardy-Weinberg P value; zero disables filtering.}
\item{pass_only}{Keep VCF records marked PASS or unfiltered.}
\item{biallelic_snps_only}{Keep biallelic single-nucleotide variants only.}
\item{backend}{VCF extraction backend.}
\item{max_variants}{Maximum retained variants.}
\item{on_excess}{Action when \code{max_variants} is exceeded.}
\item{quiet}{Suppress filtering messages.}
\item{prefix}{PLINK file prefix.}
\item{x}{Input object or path.}
\item{format}{Input format.}
\item{variants}{Variant indices, IDs, or logical selector.}
\item{...}{Arguments forwarded to a format-specific reader.}
}
\value{\code{parse_region()} returns a one-row coordinate table. Other functions return an \code{ld_data} object.}
\examples{
vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
x <- read_vcf_region(vcf, "chr1:1000000-1100000", min_maf = 0.01,
                     quiet = TRUE)
x
subset_ld_data(x, variants = 1:10, samples = 1:20)
}
