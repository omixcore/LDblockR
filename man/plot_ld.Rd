\name{plot_ld}
\alias{plot_ld}
\alias{save_ld_plot}
\title{Draw and Save Regional LD Figures}
\description{Create publication-ready inverted-triangle LD heatmaps aligned to association, gene-model, SNP-label, key-SNP, block, tag-SNP, and highlighted-variant tracks. Association rows are matched to exact LD SNP IDs or chromosome/position keys; the default layout puts the regional GWAS track at the top, optional gene/MAF tracks below it, an aligned SNP ID/position track above a key-SNP and block track immediately above the heatmap, and colored block-specific outlines plus Block labels inside the corresponding LD sub-triangles. Dashed guides connect every SNP marker directly to the actual top edge of the compact index-spaced heatmap, while custom heatmap color ramps can be supplied independently. Aligned connector lines, region labels, cell values for small regions, and a continuous LD color key are also supported. Genomic coordinates are shown on the upper GWAS track and the lower heatmap axis is suppressed. The optional MAF track is off by default; use \code{plot_ld_region()} for the reference-style integrated view with MAF and key-SNP tracks enabled.}
\usage{
plot_ld(ld, metric = "r2", gwas = NULL,
  genes = NULL, blocks = NULL, tags = NULL, lead = NULL,
  special = NULL, cutline = Inf, title = NULL, palette = "classic",
  position_scale = c("index", "physical"),
  render = c("auto", "vector", "raster"), show_values = FALSE,
  show_maf = FALSE, show_ld_colors = FALSE,
  show_snp_labels = TRUE,
  snp_label_type = c("id_position", "id", "position"),
  max_snp_labels = 120L,
  show_snp_connectors = TRUE,
  show_key_snps = TRUE,
  key_snp_colors = c(lead = "#7B2CBF", tag = "#F28E2B",
  block = "#2166AC", special = "#D73027"),
  gwas_color_by = c("significance", "pve", "ld"),
  show_pve_legend = TRUE, show_connectors = TRUE,
  show_region_label = TRUE,
  show_gene_names = TRUE,
  max_gene_lanes = 6L, point_size = 0.65, grid_color = "white",
  grid_width = 0.4,
  gene_colors = c(CDS = "#D95F02", UTR = "#7570B3",
  intron = "#80B1D3", gene = "#1B9E77", arrow = "#333333"),
  raster_width = 1400L, raster_height = 700L, draw = TRUE,
  heatmap_colors = NULL)

save_ld_plot(plot, file, width = 11, height = 8.5, dpi = 300,
  bg = "white")
}
\arguments{
\item{ld}{An \code{ld_result} object.}
\item{metric}{One or both LD metrics.}
\item{gwas}{Association file or table.}
\item{genes}{GFF3 path or annotation table.}
\item{blocks}{Optional block table. Blocks are drawn in a dedicated track immediately above the LD heatmap and labeled with colored outlines inside their corresponding LD sub-triangles; both \code{ld_blocks} output and \code{chr/start/end} tables are accepted.}
\item{tags}{Tag-SNP table.}
\item{lead}{Lead variant ID or position.}
\item{special}{Highlighted variants or named labels.}
\item{cutline}{Association minus-log10 P threshold.}
\item{title}{Figure title.}
\item{palette}{Palette name or color vector.}
\item{heatmap_colors}{Optional palette name or color vector for the LD heatmap;
when supplied, it overrides \code{palette} for the heatmap and color key.}
\item{position_scale}{Physical or equally spaced coordinates. The default is
equally spaced SNP indices shared by all tracks on the full plotting width;
this keeps the SNP/Block marker order aligned.}
\item{render}{Heatmap rendering mode.}
\item{show_values}{Print small-matrix cell values.}
\item{show_maf, show_ld_colors, show_gene_names}{Track display switches.}
\item{show_snp_labels}{Add the aligned SNP marker/label track above the Block
track and place a small marker tick for every LD SNP at the heatmap baseline.}
\item{snp_label_type}{Use SNP IDs, physical positions, or both in the label
track.}
\item{max_snp_labels}{Maximum number of labels to draw. Dense plots thin the
labels evenly; use \code{Inf} to label every SNP.}
\item{show_snp_connectors}{Draw faint dashed guides from every SNP marker
directly to the actual top edge of the LD heatmap.}
\item{show_key_snps}{Add a compact key-SNP track with lead-SNP, tag-SNP,
block-boundary, and highlighted-variant markers.}
\item{key_snp_colors}{Named colors for the \code{lead}, \code{tag},
\code{block}, and \code{special} key-SNP markers.}
\item{gwas_color_by}{Color association points by significance, PVE, or LD.}
\item{show_pve_legend}{Show the PVE point-size key when PVE is available.}
\item{show_connectors}{Draw lines from exact GWAS SNP coordinates to the
heatmap marker baseline.}
\item{show_region_label}{Show Start/Region/End text in the gene panel.}
\item{max_gene_lanes}{Maximum gene lanes.}
\item{point_size}{Association point size.}
\item{grid_color, grid_width}{Vector-cell borders.}
\item{gene_colors}{Named feature colors.}
\item{raster_width, raster_height}{Large-heatmap raster resolution.}
\item{draw}{Draw immediately.}
\item{plot}{An \code{ld_plot} object.}
\item{file}{PDF, SVG, PNG, or TIFF output path.}
\item{width, height}{Device size in inches.}
\item{dpi}{Raster resolution.}
\item{bg}{Background color.}
}
\value{\code{plot_ld()} returns an \code{ld_plot}; \code{save_ld_plot()} returns the path invisibly.}
\examples{
vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
x <- read_vcf_region(vcf, min_maf = 0.01, quiet = TRUE)
ld <- ld_compute(x, "r2")
p <- plot_ld(ld, "r2", draw = FALSE)
\dontrun{save_ld_plot(p, "LD.pdf")}
}
