#!/usr/bin/env Rscript

# Reproduce the built-in TASSEL example MLM scan and a publication-style
# regional LD figure.  Run with: Rscript scripts/reproduce_gwas_example.R out

suppressPackageStartupMessages(library(LDblockR))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "LDblockR_EarHT_example"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- example_data("tassel")
pheno <- read_tassel_phenotype(paths[["mdp_phenotype"]])
gwas <- gwas_mlm(
  genotype = paths[["mdp_genotype"]], phenotype = pheno,
  trait = "EarHT", covariates = c("location", "Q1", "Q2", "Q3"),
  replicate = "expand", min_maf = 0.05, max_missing = 0.20,
  chunk_size = 256L, verbose = TRUE
)
write.table(gwas, file.path(out_dir, "EarHT_MLM_GWAS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cutline <- -log10(0.05 / nrow(gwas))
lead0 <- gwas[which.max(gwas$logp), , drop = FALSE]
chr_levels <- unique(as.character(gwas$chr))
chr_colors <- stats::setNames(
  rep(c("#2C7FB8", "#7B3294", "#4DAF4A", "#FF7F00"),
      length.out = length(chr_levels)),
  chr_levels
)
manhattan <- plot_manhattan(
  gwas, cutline = cutline, color_by = "chromosome", chr_colors = chr_colors,
  color_cycle = TRUE, color_significant = TRUE, point_size_by = "PVE",
  highlight = lead0$id, label_top = 5L, width_mode = "length",
  title = "Ear height mixed-model GWAS", draw = FALSE
)
qq <- plot_qq(gwas, highlight = lead0$id,
              title = "Ear height mixed-model GWAS Q-Q plot", draw = FALSE)
save_gwas_plot(manhattan, file.path(out_dir, "Figure_EarHT_Manhattan.pdf"),
               width = 8.5, height = 5.4)
save_gwas_plot(manhattan, file.path(out_dir, "Figure_EarHT_Manhattan.svg"),
               width = 8.5, height = 5.4)
save_gwas_plot(manhattan, file.path(out_dir, "Figure_EarHT_Manhattan.png"),
               width = 8.5, height = 5.4, dpi = 400)
save_gwas_plot(qq, file.path(out_dir, "Figure_EarHT_QQ.pdf"),
               width = 6.8, height = 6.2)
save_gwas_plot(qq, file.path(out_dir, "Figure_EarHT_QQ.svg"),
               width = 6.8, height = 6.2)
save_gwas_plot(qq, file.path(out_dir, "Figure_EarHT_QQ.png"),
               width = 6.8, height = 6.2, dpi = 400)

region_info <- gwas_ld_region(gwas, cutline = cutline, flank = 250000,
                              min_width = 100000)
region <- region_info$region
lead <- region_info$lead_id
lead_row <- gwas[match(lead, gwas$id), , drop = FALSE]
gwas_region <- region_info$selected
write.table(gwas_region, file.path(out_dir, "EarHT_MLM_GWAS_regional.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
regional <- read_hapmap(paths[["mdp_genotype"]], region = region,
                         min_maf = 0.05, max_missing = 0.20, quiet = TRUE)
ld <- ld_compute(regional, measure = "r2", r2_method = "dosage")
blocks <- detect_ld_blocks(ld, method = "solid_spine", metric = "r2", spine_cut = 0.8)
tags <- select_tag_snps(ld, threshold = 0.8, blocks = blocks)

# The bundled maize files are genotype/phenotype examples and do not carry a
# gene annotation.  Supply a GFF3 path as the second command-line argument for
# a real annotation; otherwise the association and LD tracks remain exact.
gff <- if (length(args) >= 2L && file.exists(args[2L])) args[2L] else NULL
special <- data.frame(id = lead, pos = lead_row$pos, label = "Lead SNP",
                      stringsAsFactors = FALSE)
p <- plot_ld(ld, metric = "r2", gwas = gwas_region, genes = gff,
             blocks = blocks, tags = tags, lead = lead, special = special,
             cutline = cutline, palette = "publication",
             position_scale = "index",
             show_values = nrow(ld$data$variants) <= 40L,
             show_maf = FALSE, show_connectors = TRUE,
             gwas_color_by = "significance",
             title = "Ear height mixed-model GWAS and regional LD", draw = FALSE)
save_ld_plot(p, file.path(out_dir, "Figure_EarHT_LD.pdf"), width = 8.5, height = 7.2)
save_ld_plot(p, file.path(out_dir, "Figure_EarHT_LD.svg"), width = 8.5, height = 7.2)
save_ld_plot(p, file.path(out_dir, "Figure_EarHT_LD.png"), width = 8.5, height = 7.2, dpi = 400)

cat("Wrote", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("Lead SNP:", lead, "\n")
cat("Region:", region, "\n")
cat("Threshold:", format(cutline, digits = 6), "| threshold hits:",
    length(region_info$hit_ids), "\n")
