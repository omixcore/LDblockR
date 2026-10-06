\name{plot_ld_region}
\alias{plot_ld_region}
\title{Draw an Integrated Regional GWAS and LD View}
\description{Configure \code{plot_ld()} as a reference-style regional visualization. The returned figure aligns the regional GWAS signal, optional gene models, MAF, key SNPs, block boundaries, and the inverted-triangle LD heatmap to the exact variants stored in the \code{ld_result} object.}
\usage{
plot_ld_region(ld, metric = "r2", gwas = NULL, genes = NULL,
  blocks = NULL, tags = NULL, lead = NULL, special = NULL,
  cutline = Inf, title = "Regional GWAS and LD view",
  palette = "publication", heatmap_colors = NULL,
  position_scale = "index", render = "auto", show_maf = TRUE,
  show_ld_colors = TRUE, gwas_color_by = "ld",
  show_pve_legend = FALSE, show_snp_labels = FALSE,
  show_snp_connectors = TRUE, show_key_snps = TRUE,
  show_connectors = FALSE, show_gene_names = TRUE,
  show_region_label = TRUE, draw = TRUE, ...)
}
\arguments{
\item{ld}{An \code{ld_result} object.}
\item{metric}{One or both LD metrics.}
\item{gwas}{Optional regional GWAS file or table.}
\item{genes}{Optional GFF3/GTF file or normalized annotation table.}
\item{blocks}{Optional block table from \code{detect_ld_blocks()} or a
chromosome/start/end table.}
\item{tags}{Optional tag-SNP table from \code{select_tag_snps()}.}
\item{lead}{Lead variant ID, \code{chr:position}, or position.}
\item{special}{Optional highlighted variants or named labels.}
\item{cutline}{Regional GWAS \code{-log10(P)} threshold.}
\item{title}{Figure title.}
\item{palette}{LD heatmap palette.}
\item{heatmap_colors}{Optional independent LD heatmap palette.}
\item{position_scale}{Use equally spaced SNP indices or physical positions.}
\item{render}{Heatmap rendering mode.}
\item{show_maf}{Add the MAF track.}
\item{show_ld_colors}{Color association points by LD with the lead SNP.}
\item{gwas_color_by}{Association point color rule.}
\item{show_pve_legend}{Add the PVE size legend when available.}
\item{show_snp_labels}{Add the dense SNP label track.}
\item{show_snp_connectors}{Draw dashed SNP-to-heatmap guides.}
\item{show_key_snps}{Add the lead/tag/block-boundary key track.}
\item{show_connectors}{Add a separate vertical GWAS connector track.}
\item{show_gene_names, show_region_label}{Gene-track display switches.}
\item{draw}{Draw immediately.}
\item{...}{Additional arguments passed to \code{plot_ld()}.}
}
\value{An object of class \code{ld_plot}.}
\examples{
\dontrun{
regional <- example_data("regional")
x <- read_vcf_region(regional[["vcf"]], region = "chr1:1000000-1100000",
                     min_maf = 0.01, quiet = TRUE)
ld <- ld_compute(x, "r2")
gwas <- read_gwas(regional[["regional_gwas"]])
blocks <- read.table(regional[["fixed_blocks"]], header = TRUE,
                     sep = "\\t", stringsAsFactors = FALSE)
tags <- select_tag_snps(ld, threshold = 0.8)
p <- plot_ld_region(ld, gwas = gwas, genes = regional[["gff3"]],
                    blocks = blocks, tags = tags, lead = "snp21",
                    draw = FALSE)
save_ld_plot(p, "regional_integrated_LD.pdf", width = 8.5, height = 8)
}
}
