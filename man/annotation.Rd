\name{annotation_tracks}
\alias{read_gff3}
\alias{read_gwas}
\title{Read Gene Annotation and Association Tracks}
\description{Normalize regional GFF3/GTF features and association statistics for aligned LD plotting. A bare built-in example name such as \code{regional_gwas.tsv} is resolved from the installed package data.}
\usage{
read_gff3(file, region = NULL,
  feature_types = c("gene", "mRNA", "transcript", "exon", "CDS",
  "five_prime_UTR", "three_prime_UTR", "UTR"))

read_gwas(x, region = NULL, chr_col = NULL, pos_col = NULL,
  p_col = NULL, id_col = NULL, value_is_logp = FALSE,
  header = NULL, sep = "")
}
\arguments{
\item{file}{GFF3/GTF path, optionally compressed.}
\item{region}{Optional genomic region.}
\item{feature_types}{Feature types to retain.}
\item{x}{Association file or data frame. A bare \code{regional_gwas.tsv} name reads the bundled regional example after installation.}
\item{chr_col, pos_col, p_col, id_col}{Explicit column names or \code{NULL} for detection.}
\item{value_is_logp}{Whether the value column is already minus-log10 P.}
\item{header}{Whether an association file has a header; \code{NULL} auto-detects.}
\item{sep}{Association file separator; empty means arbitrary whitespace.}
}
\value{A normalized data frame suitable for \code{plot_ld()}.}
