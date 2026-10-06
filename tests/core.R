library(LDblockR)

vcf <- system.file("extdata", "example.vcf.gz", package = "LDblockR")
hmp <- system.file("extdata", "example.hmp.txt", package = "LDblockR")
gwas_file <- system.file("extdata", "example_gwas.tsv", package = "LDblockR")
regional_gwas_file <- system.file("extdata", "regional_gwas.tsv", package = "LDblockR")
gff_file <- system.file("extdata", "example.gff3", package = "LDblockR")
map_file <- system.file("extdata", "example_map.tsv", package = "LDblockR")
matrix_file <- system.file("extdata", "example_matrix.tsv", package = "LDblockR")

region <- "chr1:1000000-1100000"
x <- read_vcf_region(vcf, region, min_maf = 0.01, quiet = TRUE)
stopifnot(inherits(x, "ld_data"), x$n_samples == 60L, x$n_variants == 42L)
stopifnot(all(x$variants$maf >= 0.01), nrow(x$haplotypes) == 120L)

xsub <- subset_ld_data(x, variants = 1:10, samples = 1:20)
stopifnot(xsub$n_variants == 10L, xsub$n_samples == 20L)

xh <- read_hapmap(hmp, region, min_maf = 0.01, quiet = TRUE)
stopifnot(xh$n_samples == 60L, xh$n_variants == 42L)
stopifnot(all.equal(x$genotypes, xh$genotypes, check.attributes = FALSE) == TRUE)
xh_numeric_chr <- read_hapmap(hmp, "1:1000000-1100000", min_maf = 0.01, quiet = TRUE)
stopifnot(xh_numeric_chr$n_variants == 42L,
          all.equal(xh_numeric_chr$genotypes, xh$genotypes, check.attributes = FALSE) == TRUE)

mp <- read.table(map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
xm <- read_genotypes(matrix_file, format = "matrix", map = mp)
stopifnot(xm$n_samples == 60L, xm$n_variants == 42L)

# Verify native PLINK BED decoding.
plink_prefix <- file.path(tempdir(), "ldblockr_plink")
g <- matrix(c(0, 1, 2, NA, 2, 2, 1, 0, 0, 0, 1, 1), nrow = 4L)
write.table(data.frame("F", paste0("I", 1:4), 0, 0, 0, -9),
            paste0(plink_prefix, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame("1", paste0("p", 1:3), 0, 100 * 1:3, "A", "G"),
            paste0(plink_prefix, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)
con <- file(paste0(plink_prefix, ".bed"), "wb")
writeBin(as.raw(c(108, 27, 1)), con)
code <- c(`0` = 0L, `1` = 2L, `2` = 3L)
for (j in 1:3) {
  z <- ifelse(is.na(g[, j]), 1L, code[as.character(g[, j])])
  byte <- sum(bitwShiftL(as.integer(z), 2L * (0:3)))
  writeBin(as.raw(byte), con)
}
close(con)
xp <- read_plink(plink_prefix, min_maf = 0, max_missing = 1, quiet = TRUE)
stopifnot(all.equal(unname(xp$genotypes), g, check.attributes = FALSE) == TRUE)

ld <- ld_compute(x, measure = "both", r2_method = "dosage", ci = TRUE, min_n = 5L)
stopifnot(inherits(ld, "ld_result"))
stopifnot(isTRUE(all.equal(ld$r2, t(ld$r2))), isTRUE(all.equal(ld$dprime, t(ld$dprime))))
stopifnot(all(diag(ld$r2) == 1), all(diag(ld$dprime) == 1))
stopifnot(all(ld$r2[is.finite(ld$r2)] >= 0 & ld$r2[is.finite(ld$r2)] <= 1))
stopifnot(all(ld$dprime[is.finite(ld$dprime)] >= 0 & ld$dprime[is.finite(ld$dprime)] <= 1))
manual_r2 <- suppressWarnings(cor(x$genotypes, use = "pairwise.complete.obs")^2)
stopifnot(max(abs(ld$r2 - manual_r2), na.rm = TRUE) < 1e-12)
ld_auto <- ld_compute(x, measure = "both")
manual_haplotype_r2 <- suppressWarnings(cor(x$haplotypes, use = "pairwise.complete.obs")^2)
stopifnot(ld_auto$r2_method == "haplotype")
stopifnot(max(abs(ld_auto$r2 - manual_haplotype_r2), na.rm = TRUE) < 1e-12)

gabriel <- detect_ld_blocks(ld, "gabriel")
strong <- detect_ld_blocks(ld, "strong", metric = "r2", strong_cut = 0.7,
                           strong_fraction = 0.7)
spine <- detect_ld_blocks(ld, "solid_spine", metric = "r2", spine_cut = 0.65)
fg <- detect_ld_blocks(ld, "four_gamete")
fixed_input <- read.table(system.file("extdata", "example_fixed_blocks.tsv", package = "LDblockR"),
                          header = TRUE)
fixed <- detect_ld_blocks(ld, "fixed", fixed = fixed_input)
stopifnot(inherits(gabriel, "ld_blocks"), nrow(strong) > 0, nrow(spine) > 0,
          inherits(fg, "ld_blocks"), nrow(fixed) == 3L)

tags <- select_tag_snps(ld, threshold = 0.8, blocks = strong)
stopifnot(nrow(tags) > 0L, all(tags$tag %in% x$variants$id))
global_tags <- select_tag_snps(ld, threshold = 0.8)
covered <- unique(unlist(strsplit(global_tags$covered, "|", fixed = TRUE)))
stopifnot(setequal(covered, x$variants$id))

neighbors <- ld_neighbors(ld, "snp21", threshold = 0.5)
stopifnot("snp21" %in% neighbors$id)
haps <- haplotype_frequencies(x, 1:5)
stopifnot(abs(sum(haps$frequency) - 1) < 0.1, all(haps$frequency >= 0.01))
block_haps <- haplotype_frequencies(x, "snp1|snp2|snp3")
stopifnot(identical(attr(block_haps, "variants"), c("snp1", "snp2", "snp3")))
empty_variant_message <- tryCatch(haplotype_frequencies(x, character()),
                                  error = function(e) conditionMessage(e))
stopifnot(is.character(empty_variant_message),
          grepl("No variants were supplied", empty_variant_message, fixed = TRUE))
decay <- ld_decay(ld, "r2", n_bins = 12L)
stopifnot(nrow(decay) > 1L, all(decay$n_pairs > 0L))

groups <- setNames(rep(c("A", "B"), each = 30L), x$samples)
group_ld <- ld_by_group(x, groups, measure = "r2")
delta <- ld_compare(group_ld, "r2", "A")
stopifnot(length(group_ld) == 2L, identical(dim(delta[[1L]]), c(42L, 42L)))

gwas <- read_gwas(gwas_file, region)
regional_paths <- example_data("regional")
regional_gwas <- read_gwas("regional_gwas.tsv", region)
genes <- read_gff3(gff_file, region)
special <- read.table(system.file("extdata", "example_special.tsv", package = "LDblockR"), header = TRUE)
stopifnot(file.exists(regional_gwas_file), file.exists(regional_paths[["regional_gwas"]]),
          nrow(gwas) == 42L, nrow(regional_gwas) == 42L, nrow(genes) > 3L)
gwas_region <- gwas_ld_region(gwas, cutline = 5, flank = 5000, min_width = 10000)
stopifnot(inherits(gwas_region, "gwas_region"), nchar(gwas_region$region) > 0,
          all(gwas_region$selected$chr == gwas_region$chr),
          all(gwas_region$selected$pos >= gwas_region$start),
          all(gwas_region$selected$pos <= gwas_region$end))

p <- plot_ld(ld, "both", gwas = gwas, genes = genes, blocks = strong,
             tags = tags, special = special, cutline = 5, draw = FALSE)
plot_files <- file.path(tempdir(), c("ldplot.pdf", "ldplot.svg", "ldplot.png", "ldplot.tiff"))
for (f in plot_files) save_ld_plot(p, f, width = 7, height = 7, dpi = 90)
stopifnot(all(file.exists(plot_files)), all(file.info(plot_files)$size > 1000L))
fixed_blocks <- read.table(regional_paths[["fixed_blocks"]], header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE)
fixed_tags <- select_tag_snps(ld, threshold = 0.8, blocks = fixed_blocks)
fixed_plot <- plot_ld(ld, "r2", gwas = regional_gwas, blocks = fixed_blocks,
                      draw = FALSE)
fixed_plot_file <- file.path(tempdir(), "ldplot_fixed_blocks.png")
save_ld_plot(fixed_plot, fixed_plot_file, width = 7, height = 5.5, dpi = 90)
block_coord <- LDblockR:::.plot_coordinates(x$variants, "index")
block_spans <- LDblockR:::.plot_block_spans(fixed_blocks, block_coord)
snp_labels <- LDblockR:::.snp_plot_labels(x$variants, "id_position")
stopifnot(file.exists(fixed_plot_file), file.info(fixed_plot_file)$size > 1000L,
          nrow(block_spans) == 3L,
          identical(block_spans$label, paste0("Block", 1:3)),
          all(block_spans$right > block_spans$left),
          length(snp_labels) == nrow(x$variants),
          all(grepl("\n", snp_labels, fixed = TRUE)),
          nrow(fixed_tags) > 0L,
          all(grepl("^Block[123]$|^Unblocked_", fixed_tags$block)),
          identical(LDblockR:::.snp_plot_labels(x$variants, "id"),
                    as.character(x$variants$id)),
          identical(LDblockR:::.snp_plot_labels(x$variants, "position"),
                    format(x$variants$pos, scientific = FALSE, trim = TRUE,
                           big.mark = ",")),
          isTRUE(fixed_plot$show_snp_labels),
          identical(fixed_plot$snp_label_type, "id_position"),
          isTRUE(fixed_plot$show_snp_connectors))
snp_id_plot <- plot_ld(ld, "r2", show_snp_labels = TRUE,
                       snp_label_type = "id", max_snp_labels = Inf,
                       draw = FALSE)
snp_hidden_plot <- plot_ld(ld, "r2", show_snp_labels = FALSE,
                           show_snp_connectors = FALSE, draw = FALSE)
custom_color_plot <- plot_ld(ld, "r2", palette = "red",
                             heatmap_colors = c("white", "skyblue", "navy"),
                             draw = FALSE)
integrated_plot <- plot_ld_region(
  ld, metric = "r2", gwas = gwas, genes = genes, blocks = strong,
  tags = tags, lead = "snp21", special = special, cutline = 5,
  draw = FALSE
)
stopifnot(identical(snp_id_plot$snp_label_type, "id"),
          isTRUE(is.infinite(snp_id_plot$max_snp_labels)),
          !isTRUE(snp_hidden_plot$show_snp_labels),
          !isTRUE(snp_hidden_plot$show_snp_connectors),
          is.character(custom_color_plot$heatmap_colors),
          length(custom_color_plot$heatmap_colors) == 101L,
          identical(custom_color_plot$palette, custom_color_plot$heatmap_colors),
          custom_color_plot$heatmap_colors[1L] != fixed_plot$palette[1L],
          inherits(integrated_plot, "ld_plot"),
          isTRUE(integrated_plot$show_maf),
          isTRUE(integrated_plot$show_key_snps),
          identical(integrated_plot$gwas_color_by, "ld"),
          identical(integrated_plot$show_connectors, FALSE),
          identical(names(integrated_plot$key_snp_colors),
                    c("lead", "tag", "block", "special")))

exported <- export_ld(ld, file.path(tempdir(), "ldstats"), strong, tags, compatible = TRUE)
stopifnot(length(exported) == 8L, all(file.exists(exported)), all(file.info(exported)$size > 0L))
triangle_lines <- readLines(gzfile(exported[["triangle_compatible"]]))
stopifnot(length(triangle_lines) == 42L,
          length(strsplit(triangle_lines[2L], "\t", fixed = TRUE)[[1L]]) == 41L)
site_first <- readLines(gzfile(exported[["site_compatible"]]), n = 1L)
stopifnot(startsWith(site_first, "chr1\t"), !grepl("pos", site_first, fixed = TRUE))

# Force and verify raster rendering on a larger matrix.
set.seed(1)
large_g <- matrix(rbinom(24L * 260L, 2L, 0.3), nrow = 24L)
large_map <- data.frame(chr = "1", pos = seq_len(260L) * 100L, id = paste0("v", 1:260))
large_x <- as_ld_data(large_g, large_map)
large_ld <- ld_compute(large_x, "r2")
large_plot <- plot_ld(large_ld, "r2", render = "auto", show_maf = FALSE, draw = FALSE)
large_png <- file.path(tempdir(), "large_raster.png")
save_ld_plot(large_plot, large_png, width = 7, height = 4, dpi = 80)
stopifnot(file.exists(large_png), file.info(large_png)$size > 1000L)

# Full one-call workflow.
prefix <- file.path(tempdir(), "workflow", "demo")
result <- ldblockr(vcf, prefix, region = region, min_maf = 0.01,
                   measure = "both", block_method = "strong",
                   block_metric = "r2", gwas = gwas_file, gff = gff_file,
                   output_formats = "pdf",
                   block_args = list(strong_cut = 0.7, strong_fraction = 0.7))
stopifnot(inherits(result$ld, "ld_result"), all(file.exists(result$files)))
stopifnot(all(c("sites", "pairs", "plot_pdf") %in% names(result$files)))

matrix_prefix <- file.path(tempdir(), "workflow", "matrix_demo")
matrix_result <- ldblockr(matrix_file, matrix_prefix, region = region,
                          format = "matrix", read_args = list(map = mp),
                          min_maf = 0.01, measure = "r2", block_method = "none",
                          output_formats = "pdf")
stopifnot(matrix_result$data$n_variants == 42L, all(file.exists(matrix_result$files)))

# Fast mixed-model GWAS smoke test and PVE contract.
set.seed(17)
n_g <- 36L
n_m <- 48L
g_gwas <- matrix(rbinom(n_g * n_m, 2L, 0.35), nrow = n_g,
                 dimnames = list(paste0("taxon", seq_len(n_g)), paste0("m", seq_len(n_m))))
g_map <- data.frame(chr = "1", pos = seq_len(n_m) * 1000, id = colnames(g_gwas))
p_gwas <- data.frame(Taxa = rownames(g_gwas), location = rep(c("A", "B"), each = n_g / 2),
                      trait = as.numeric(g_gwas[, 3] + rnorm(n_g)), Q1 = runif(n_g),
                      stringsAsFactors = FALSE)
gwas_fit <- gwas_mlm(g_gwas, p_gwas, trait = "trait", map = g_map,
                     covariates = c("location", "Q1"), min_maf = 0.05,
                     max_missing = 0.1, chunk_size = 12L, verbose = FALSE)
stopifnot(inherits(gwas_fit, "gwas_mlm_result"), nrow(gwas_fit) > 0L,
          all(c("p", "PVE", "model_PVE") %in% names(gwas_fit)),
          all(gwas_fit$PVE >= 0 & gwas_fit$PVE <= 1, na.rm = TRUE),
          is.finite(attr(gwas_fit, "model")$model_PVE))
# Regression checks for the publication workflow: subsetting/head() must not
# recurse through the custom result printer, and metric omission must default
# to the r2 heatmap used in the manuscript figure.
preview <- capture.output(head(gwas_fit[order(gwas_fit$p),
                                         c("id", "chr", "pos", "p", "PVE")]))
stopifnot(length(preview) > 0L,
          any(grepl("gwas_mlm_result", preview, fixed = TRUE)))
# Manhattan and Q-Q diagnostics use the same normalized GWAS result and can be
# saved in every publication format without external plotting dependencies.
manhattan <- plot_manhattan(gwas_fit, color_by = "significance", label_top = 2L,
                            draw = FALSE)
qq <- plot_qq(gwas_fit, draw = FALSE)
stopifnot(inherits(manhattan, "manhattan_plot"), inherits(qq, "qq_plot"),
          length(manhattan$data$x) == nrow(gwas_fit),
          length(qq$expected) <= nrow(gwas_fit), is.finite(qq$lambda),
          manhattan$cutline > 0)
# CMplot-style controls: named per-chromosome colors, optional color cycling,
# multiple thresholds, custom labels, and significance-size emphasis.
multi_gwas <- as.data.frame(gwas_fit, stringsAsFactors = FALSE)
multi_gwas$chr <- rep(c("1", "2", "3"), length.out = nrow(multi_gwas))
multi_gwas$pos <- seq_len(nrow(multi_gwas)) * 1000
multi_named <- plot_manhattan(
  multi_gwas, threshold = c(2, 4), suggestive = 1,
  chr_colors = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3"),
  color_cycle = FALSE, chr_labels = c("Chr I", "Chr II", "Chr III"),
  color_significant = TRUE, pch = 17, signal_cex = 1.4, alpha = 0.7,
  draw = FALSE)
multi_cycle <- plot_manhattan(
  multi_gwas, threshold = 3, chr_colors = c("#111111", "#BBBBBB"),
  color_cycle = TRUE, draw = FALSE)
m_width <- plot_manhattan(
  multi_gwas, chr_lengths = c("1" = 100000, "2" = 150000, "3" = 200000),
  width_mode = "length", chr_colors = c("1" = "#1B9E77",
                                         "2" = "#D95F02", "3" = "#7570B3"),
  color_cycle = FALSE, draw = FALSE)
gap_fixed <- plot_manhattan(
  multi_gwas, chr_lengths = c("1" = 100000, "2" = 150000, "3" = 200000),
  width_mode = "length", gap_bp = 5000, gap_mode = "fixed", draw = FALSE)
gap_relative <- plot_manhattan(
  multi_gwas, chr_lengths = c("1" = 100000, "2" = 150000, "3" = 200000),
  width_mode = "length", gap_ratio = 0.05, gap_mode = "relative", draw = FALSE)
stopifnot(identical(multi_named$chr_labels, c("Chr I", "Chr II", "Chr III")),
          length(multi_named$threshold) == 2L,
          length(multi_named$point_colors) == nrow(multi_gwas),
          length(unique(multi_cycle$chr)) == 3L,
          length(multi_cycle$point_colors) == nrow(multi_gwas),
          identical(as.numeric(m_width$chr_lengths), c(100000, 150000, 200000)),
          m_width$chr_spans[1L] < m_width$chr_spans[2L],
          m_width$chr_spans[2L] < m_width$chr_spans[3L],
          min(m_width$chr_mid - m_width$chr_spans / 2) <= min(m_width$data$x),
          max(m_width$chr_mid + m_width$chr_spans / 2) > max(m_width$data$x),
          identical(gap_fixed$gap_mode, "fixed"),
          length(gap_fixed$gap_values) == 2L,
          all(gap_fixed$gap_values == 5000),
          identical(gap_relative$gap_mode, "relative"),
          length(gap_relative$gap_values) == 2L,
          !identical(gap_relative$gap_values[1L], gap_relative$gap_values[2L]))
diag_files <- file.path(tempdir(), c("manhattan.pdf", "manhattan.png",
                                     "qq.svg", "qq.tiff"))
save_gwas_plot(manhattan, diag_files[1L], width = 6, height = 4.5)
save_gwas_plot(manhattan, diag_files[2L], width = 6, height = 4.5, dpi = 90)
save_gwas_plot(qq, diag_files[3L], width = 5.5, height = 5.5)
save_gwas_plot(qq, diag_files[4L], width = 5.5, height = 5.5, dpi = 90)
stopifnot(all(file.exists(diag_files)), all(file.info(diag_files)$size > 1000L))
default_metric_plot <- plot_ld(ld, gwas = gwas, draw = FALSE)
stopifnot(inherits(default_metric_plot, "ld_plot"),
          identical(default_metric_plot$metrics, "r2"))

tassel_paths <- example_data("tassel")
stopifnot(file.exists(tassel_paths[["mdp_genotype"]]),
          basename(tassel_paths[["mdp_genotype"]]) == "mdp_genotype.hmp.txt",
          file.info(tassel_paths[["mdp_genotype"]])$size > 1000000L)
tassel_pheno <- read_tassel_phenotype(tassel_paths[["mdp_phenotype"]])
stopifnot(nrow(tassel_pheno) > 500L, anyDuplicated(tassel_pheno$Taxa) > 0L)
tassel_k <- read_tassel_kinship(tassel_paths[["mdp_kinship"]])
stopifnot(nrow(tassel_k) == ncol(tassel_k), !is.null(rownames(tassel_k)))
gwas_norm <- read_gwas(data.frame(chr = "1", pos = 10, p = 0.01,
                                  id = "m1", PVE = 0.12))
stopifnot("PVE" %in% names(gwas_norm), identical(gwas_norm$PVE, 0.12))
