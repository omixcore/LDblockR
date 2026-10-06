#!/usr/bin/env Rscript

# Reproduce the reference-style integrated regional view from the bundled
# synthetic regional fixture. Run with:
#   Rscript scripts/reproduce_integrated_region.R integrated_region_results

suppressPackageStartupMessages(library(LDblockR))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "LDblockR_integrated_region"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- example_data("regional")
region <- "chr1:1000000-1100000"
x <- read_vcf_region(paths[["vcf"]], region = region,
                     min_maf = 0.01, max_missing = 0.25, quiet = TRUE)
ld <- ld_compute(x, measure = "r2", r2_method = "dosage")
gwas <- read_gwas(paths[["regional_gwas"]], region = region)
blocks <- read.table(paths[["fixed_blocks"]], header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)
tags <- select_tag_snps(ld, threshold = 0.8, blocks = blocks)
special <- read.table(paths[["special"]], header = TRUE,
                      sep = "\t", stringsAsFactors = FALSE)
lead <- gwas$id[which.max(gwas$logp)]
cutline <- 5

p <- plot_ld_region(
  ld, gwas = gwas, genes = paths[["gff3"]], blocks = blocks,
  tags = tags, lead = lead, special = special, cutline = cutline,
  heatmap_colors = c("#FFFDF2", "#FDBB55", "#B40426"),
  show_values = TRUE, draw = FALSE
)

save_ld_plot(p, file.path(out_dir, "Figure_integrated_regional_LD.pdf"),
             width = 8.5, height = 8)
save_ld_plot(p, file.path(out_dir, "Figure_integrated_regional_LD.svg"),
             width = 8.5, height = 8)
save_ld_plot(p, file.path(out_dir, "Figure_integrated_regional_LD.png"),
             width = 8.5, height = 8, dpi = 400)
save_ld_plot(p, file.path(out_dir, "Figure_integrated_regional_LD.tiff"),
             width = 8.5, height = 8, dpi = 400)

write.table(gwas, file.path(out_dir, "regional_gwas.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(blocks, file.path(out_dir, "blocks.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(tags, file.path(out_dir, "tags.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Wrote", normalizePath(out_dir, mustWork = FALSE), "\n")
