.write_table_connection <- function(x, file, col.names = TRUE) {
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "wt") else base::file(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = col.names, na = "NA")
  invisible(file)
}

.write_pairwise <- function(ld, file, compatible = FALSE) {
  con <- gzfile(file, "wt")
  on.exit(close(con), add = TRUE)
  v <- ld$data$variants
  wrote_header <- FALSE
  for (i in seq_len(nrow(v) - 1L)) {
    js <- (i + 1L):nrow(v)
    available <- rep(FALSE, length(js))
    for (metric in c("r2", "dprime")) {
      if (!is.null(ld[[metric]])) available <- available | is.finite(ld[[metric]][i, js])
    }
    js <- js[available]
    if (!length(js)) next
    if (compatible) {
      z <- data.frame(
        chr = v$chr[i], pos1 = v$pos[i], pos2 = v$pos[js],
        r2 = if (is.null(ld$r2)) NA_real_ else ld$r2[i, js],
        dprime = if (is.null(ld$dprime)) NA_real_ else ld$dprime[i, js],
        stringsAsFactors = FALSE
      )
    } else {
      z <- data.frame(
        chr = v$chr[i], pos1 = v$pos[i], pos2 = v$pos[js],
        id1 = v$id[i], id2 = v$id[js], distance = v$pos[js] - v$pos[i],
        n = ld$n[i, js], r2 = if (is.null(ld$r2)) NA_real_ else ld$r2[i, js],
        dprime = if (is.null(ld$dprime)) NA_real_ else ld$dprime[i, js],
        D = if (is.null(ld$D)) NA_real_ else ld$D[i, js],
        ci_lower = if (is.null(ld$ci_lower)) NA_real_ else ld$ci_lower[i, js],
        ci_upper = if (is.null(ld$ci_upper)) NA_real_ else ld$ci_upper[i, js],
        stringsAsFactors = FALSE
      )
    }
    utils::write.table(z, con, sep = "\t", quote = FALSE, row.names = FALSE,
                       col.names = !wrote_header, na = "NA")
    wrote_header <- TRUE
  }
  if (!wrote_header) {
    z <- if (compatible) data.frame(chr = character(), pos1 = numeric(), pos2 = numeric(),
                                    r2 = numeric(), dprime = numeric()) else
      data.frame(chr = character(), pos1 = numeric(), pos2 = numeric(), id1 = character(),
                 id2 = character(), distance = numeric(), n = integer(), r2 = numeric(),
                 dprime = numeric(), D = numeric(), ci_lower = numeric(), ci_upper = numeric())
    utils::write.table(z, con, sep = "\t", quote = FALSE, row.names = FALSE,
                       col.names = TRUE, na = "NA")
  }
  invisible(file)
}

.write_compat_triangle <- function(mat, file, label) {
  con <- gzfile(file, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("#Region PairWise ", label), con)
  if (nrow(mat) >= 2L) {
    for (i in seq_len(nrow(mat) - 1L)) {
      values <- mat[i, (i + 1L):ncol(mat)]
      text <- ifelse(is.finite(values), formatC(values, format = "f", digits = 3), "NA")
      writeLines(paste(text, collapse = "\t"), con)
    }
  }
  invisible(file)
}

.write_compat_sites <- function(variants, file) {
  con <- gzfile(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(variants[, c("chr", "pos"), drop = FALSE], con, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA")
  invisible(file)
}

.write_compat_blocks <- function(blocks, file) {
  con <- gzfile(file, "wt")
  on.exit(close(con), add = TRUE)
  writeLines("#chr\tStart\tEnd\tSNPNumber\tTagSNPList", con)
  if (!is.null(blocks) && nrow(blocks)) {
    z <- data.frame(chr = blocks$chr, Start = blocks$start, End = blocks$end,
                    SNPNumber = blocks$n_snps, TagSNPList = blocks$snps,
                    stringsAsFactors = FALSE)
    utils::write.table(z, con, sep = "\t", quote = FALSE, row.names = FALSE,
                       col.names = FALSE, na = "NA")
  }
  invisible(file)
}

#' Export LD statistics and LDBlockShow-compatible tables
#'
#' @param ld An `ld_result` object.
#' @param prefix Output prefix.
#' @param blocks Optional blocks.
#' @param tags Optional tag SNPs.
#' @param compatible Also create `.site.gz`, `.blocks.gz`, and `.TriangleV.gz` files.
#' @return Named output paths.
#' @export
export_ld <- function(ld, prefix, blocks = NULL, tags = NULL, compatible = TRUE) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  v <- ld$data$variants
  sites <- paste0(prefix, ".sites.tsv.gz")
  pairs <- paste0(prefix, ".pairs.tsv.gz")
  .write_table_connection(v, sites)
  .write_pairwise(ld, pairs)
  paths <- c(sites = sites, pairs = pairs)
  if (!is.null(blocks)) {
    block_file <- paste0(prefix, ".blocks.tsv.gz")
    .write_table_connection(as.data.frame(blocks), block_file)
    paths <- c(paths, blocks = block_file)
  }
  if (!is.null(tags)) {
    tag_file <- paste0(prefix, ".tags.tsv")
    .write_table_connection(as.data.frame(tags), tag_file)
    paths <- c(paths, tags = tag_file)
  }
  if (compatible) {
    csite <- paste0(prefix, ".site.gz")
    cpair <- paste0(prefix, ".TriangleV.gz")
    cblock <- paste0(prefix, ".blocks.gz")
    .write_compat_sites(v, csite)
    primary <- if (identical(ld$measure, "r2") || is.null(ld$dprime)) "r2" else "dprime"
    primary_label <- if (primary == "r2") "R^2" else "D'"
    .write_compat_triangle(ld[[primary]], cpair, primary_label)
    .write_compat_blocks(blocks, cblock)
    paths <- c(paths, site_compatible = csite, triangle_compatible = cpair,
               blocks_compatible = cblock)
    secondary <- if (primary == "r2") "dprime" else "r2"
    if (identical(ld$measure, "both") && !is.null(ld[[secondary]])) {
      cb <- paste0(prefix, ".TriangleB.gz")
      .write_compat_triangle(ld[[secondary]], cb,
                             if (secondary == "r2") "R^2" else "D'")
      paths <- c(paths, triangle_secondary_compatible = cb)
    }
  }
  paths <- stats::setNames(normalizePath(unname(paths), mustWork = FALSE), names(paths))
  paths
}

.infer_input_format <- function(input, format) {
  if (format != "auto") return(format)
  if (inherits(input, "ld_data") || is.matrix(input) || is.data.frame(input)) return("matrix")
  low <- tolower(as.character(input))
  if (grepl("\\.vcf(?:\\.gz|\\.bgz)?$", low)) "vcf"
  else if (grepl("(?:hmp|hapmap)(?:\\.txt)?(?:\\.gz)?$", low)) "hapmap"
  else if (file.exists(paste0(input, ".bed")) || file.exists(paste0(input, ".ped"))) "plink"
  else "matrix"
}

#' Run the complete regional LD workflow
#'
#' @param input Genotype input accepted by `read_genotypes`.
#' @param output_prefix Output path prefix.
#' @param region Optional genomic region.
#' @param format Input format.
#' @param samples Optional sample subset.
#' @param min_maf,max_missing,max_het,min_hwe_p Variant filters.
#' @param measure LD statistic selection.
#' @param r2_method Automatic, dosage, or phased-haplotype r-squared.
#' @param max_distance Maximum calculated pair distance.
#' @param min_n Minimum pairwise sample count.
#' @param block_method Block algorithm.
#' @param block_metric Block metric for applicable algorithms.
#' @param tag_threshold Tag-SNP threshold.
#' @param gwas Optional GWAS file/table.
#' @param gff Optional GFF3 file/table.
#' @param special Optional highlighted variants.
#' @param output_formats Any of pdf, svg, png, or tiff.
#' @param plot_metric Plot metric or both.
#' @param width,height,dpi Figure settings.
#' @param compatible Create LDBlockShow-compatible tables.
#' @param read_args,ld_args,block_args,plot_args Additional named argument lists.
#' @return A list containing data, LD, blocks, tags, plot, and file paths.
#' @export
ldblockr <- function(input, output_prefix, region = NULL,
                     format = c("auto", "vcf", "hapmap", "plink", "matrix"),
                     samples = NULL, min_maf = 0.05, max_missing = 0.25,
                     max_het = 1, min_hwe_p = 0,
                     measure = c("both", "r2", "dprime"),
                     r2_method = c("auto", "dosage", "haplotype"),
                     max_distance = Inf, min_n = 5L,
                     block_method = c("gabriel", "solid_spine", "strong", "four_gamete", "fixed", "none"),
                     block_metric = c("dprime", "r2"), tag_threshold = 0.8,
                     gwas = NULL, gff = NULL, special = NULL,
                     output_formats = c("pdf", "svg"), plot_metric = measure,
                     width = 11, height = 8.5, dpi = 300, compatible = TRUE,
                     read_args = list(), ld_args = list(), block_args = list(),
                     plot_args = list()) {
  format <- match.arg(format)
  measure <- match.arg(measure)
  r2_method <- match.arg(r2_method)
  block_method <- match.arg(block_method)
  block_metric <- match.arg(block_metric)
  actual_format <- .infer_input_format(input, format)
  if (inherits(input, "ld_data")) {
    data <- input
    if (!is.null(samples)) data <- subset_ld_data(data, samples = .read_sample_ids(samples))
    if (!is.null(region)) {
      reg <- parse_region(region)
      keep <- .in_region(data$variants$chr, data$variants$pos, reg)
      if (!any(keep)) .stopf("No variants occur in the requested region.")
      data <- subset_ld_data(data, variants = keep)
    }
    data <- filter_variants(data, min_maf = min_maf, max_missing = max_missing,
                            max_het = max_het, min_hwe_p = min_hwe_p)
  } else if (is.matrix(input) || is.data.frame(input)) {
    data <- do.call(as_ld_data, utils::modifyList(list(genotypes = input), read_args))
    data <- filter_variants(data, min_maf = min_maf, max_missing = max_missing,
                            max_het = max_het, min_hwe_p = min_hwe_p)
  } else {
    reader <- switch(actual_format, vcf = read_vcf_region, hapmap = read_hapmap,
                     plink = read_plink, matrix = read_genotypes)
    common <- list(region = region, samples = samples, min_maf = min_maf,
                   max_missing = max_missing, max_het = max_het,
                   min_hwe_p = min_hwe_p)
    if (actual_format == "plink") common$prefix <- input else common$file <- input
    if (actual_format == "matrix") {
      common <- utils::modifyList(list(x = input, format = "matrix"), read_args)
    } else {
      common <- utils::modifyList(common, read_args)
    }
    data <- do.call(reader, common)
    if (actual_format == "matrix") {
      if (!is.null(samples)) data <- subset_ld_data(data, samples = .read_sample_ids(samples))
      if (!is.null(region)) {
        reg <- parse_region(region)
        keep <- .in_region(data$variants$chr, data$variants$pos, reg)
        if (!any(keep)) .stopf("No matrix variants occur in the requested region.")
        data <- subset_ld_data(data, variants = keep)
      }
      data <- filter_variants(data, min_maf = min_maf, max_missing = max_missing,
                              max_het = max_het, min_hwe_p = min_hwe_p)
    }
  }
  effective_measure <- measure
  plot_needed <- if (length(plot_metric) == 1L && plot_metric == "both") c("r2", "dprime") else plot_metric
  block_needed <- if (block_method %in% c("strong", "solid_spine")) block_metric else if (block_method == "gabriel") "dprime" else character()
  available_needed <- unique(c(if (measure == "both") c("r2", "dprime") else measure,
                               plot_needed, block_needed))
  if (all(c("r2", "dprime") %in% available_needed)) effective_measure <- "both"
  compute_args <- utils::modifyList(list(x = data, measure = effective_measure, r2_method = r2_method,
                                         max_distance = max_distance, min_n = min_n,
                                         ci = block_method == "gabriel"), ld_args)
  ld <- do.call(ld_compute, compute_args)
  blocks <- do.call(detect_ld_blocks,
                    utils::modifyList(list(ld = ld, method = block_method, metric = block_metric), block_args))
  tag_metric <- if (!is.null(ld$r2)) "r2" else "dprime"
  tags <- select_tag_snps(ld, threshold = tag_threshold, metric = tag_metric,
                          blocks = blocks)
  region_string <- paste0(data$variants$chr[1L], ":", min(data$variants$pos), "-", max(data$variants$pos))
  if (is.character(gwas) && length(gwas) == 1L) gwas_data <- read_gwas(gwas, region_string) else gwas_data <- gwas
  if (is.character(gff) && length(gff) == 1L) gene_data <- read_gff3(gff, region_string) else gene_data <- gff
  if (length(plot_metric) == 1L && plot_metric == "both") plot_metric <- c("r2", "dprime")
  p <- do.call(plot_ld, utils::modifyList(list(ld = ld, metric = plot_metric, gwas = gwas_data,
                                               genes = gene_data, blocks = blocks, tags = tags,
                                               special = special, draw = FALSE), plot_args))
  paths <- export_ld(ld, output_prefix, blocks, tags, compatible)
  output_formats <- unique(tolower(output_formats))
  bad <- setdiff(output_formats, c("pdf", "svg", "png", "tif", "tiff"))
  if (length(bad)) .stopf("Unsupported output formats: %s", paste(bad, collapse = ", "))
  for (ext in output_formats) {
    f <- paste0(output_prefix, ".", ext)
    save_ld_plot(p, f, width = width, height = height, dpi = dpi)
    paths <- c(paths, setNames(normalizePath(f, mustWork = FALSE), paste0("plot_", ext)))
  }
  out <- list(data = data, ld = ld, blocks = blocks, tags = tags,
              plot = p, files = paths, call = match.call())
  class(out) <- "ldblockr_result"
  out
}

#' Run the complete workflow over multiple regions
#'
#' @param input Genotype input.
#' @param regions Character regions or a chr/start/end data frame.
#' @param output_dir Output directory.
#' @param continue_on_error Continue and retain error objects.
#' @param ... Arguments passed to `ldblockr`.
#' @return A named list of workflow results.
#' @export
ld_batch <- function(input, regions, output_dir, continue_on_error = TRUE, ...) {
  if (is.data.frame(regions)) {
    chr_col <- .match_column(regions, c("chr", "chrom", "chromosome"), TRUE, "chromosome")
    start_col <- .match_column(regions, c("start", "bp1", "from"), TRUE, "start")
    end_col <- .match_column(regions, c("end", "bp2", "to"), TRUE, "end")
    regions <- paste0(regions[[chr_col]], ":", regions[[start_col]], "-", regions[[end_col]])
  }
  regions <- as.character(regions)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out <- setNames(vector("list", length(regions)), regions)
  for (i in seq_along(regions)) {
    safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", regions[i])
    prefix <- file.path(output_dir, safe_name)
    ans <- try(ldblockr(input, prefix, region = regions[i], ...), silent = TRUE)
    if (inherits(ans, "try-error") && !continue_on_error) stop(ans)
    out[[i]] <- ans
  }
  class(out) <- c("ld_batch_result", "list")
  out
}

.parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (!startsWith(token, "-")) .stopf("Unexpected command-line argument: %s", token)
    key <- sub("^-+", "", token)
    if (i == length(args) || startsWith(args[i + 1L], "-")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[i + 1L]
      i <- i + 2L
    }
  }
  out
}

.cli_first <- function(x, keys, default = NULL) {
  for (key in keys) if (!is.null(x[[key]])) return(x[[key]])
  default
}

.cli_help <- function() {
  cat(paste0(
    "LDblockR command line\n\n",
    "Usage:\n",
    "  LDblockR --vcf input.vcf.gz --out result --region chr1:1000-2000 [options]\n\n",
    "Main options:\n",
    "  --vcf FILE | --plink PREFIX | --hapmap FILE\n",
    "  --out PREFIX             Output prefix\n",
    "  --region REGION          chr:start-end\n",
    "  --samples FILE           One sample ID per line\n",
    "  --measure both|r2|dprime [both]\n",
    "  --block-method gabriel|solid_spine|strong|four_gamete|fixed|none [gabriel]\n",
    "  --gwas FILE              chr position P-value table\n",
    "  --gff FILE               GFF3/GTF annotation\n",
    "  --maf FLOAT              Minimum MAF [0.05]\n",
    "  --missing FLOAT          Maximum missing fraction [0.25]\n",
    "  --het FLOAT              Maximum heterozygous fraction [1]\n",
    "  --tag-r2 FLOAT           Tag-SNP threshold [0.8]\n",
    "  --formats pdf,svg,png    Plot formats [pdf,svg]\n\n",
    "LDBlockShow aliases are accepted: -InVCF, -OutPut, -Region, -SubPop,\n",
    "-SeleVar, -BlockType, -InGWAS, -InGFF, -MAF, -Miss, -Het,\n",
    "-TagSNPCut, -OutPdf, and -OutPng.\n"
  ))
}

#' Command-line entry point
#'
#' @param args Command-line arguments, defaulting to `commandArgs(TRUE)`.
#' @return The workflow result invisibly.
#' @export
ldblockr_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .parse_cli_args(args)
  if (!length(opt) || isTRUE(.cli_first(opt, c("help", "h"), FALSE))) {
    .cli_help()
    return(invisible(NULL))
  }
  vcf <- .cli_first(opt, c("vcf", "InVCF", "i"))
  plink <- .cli_first(opt, c("plink", "InPlink"))
  hapmap <- .cli_first(opt, c("hapmap"))
  present <- !vapply(list(vcf, plink, hapmap), is.null, logical(1L))
  if (sum(present) != 1L) .stopf("Specify exactly one of --vcf, --plink, or --hapmap.")
  input <- list(vcf, plink, hapmap)[[which(present)]]
  format <- c("vcf", "plink", "hapmap")[which(present)]
  out <- .cli_first(opt, c("out", "OutPut", "o"))
  region <- .cli_first(opt, c("region", "Region", "r"))
  if (is.null(out) || is.null(region)) .stopf("--out and --region are required.")
  sele <- as.integer(.cli_first(opt, c("SeleVar"), 0L))
  measure <- .cli_first(opt, c("measure"), if (sele == 1L) "dprime" else if (sele == 2L) "r2" else "both")
  bt <- as.integer(.cli_first(opt, c("BlockType"), 0L))
  block_method <- .cli_first(opt, c("block-method"),
                             if (bt == 2L) "solid_spine" else if (bt == 3L) "strong" else if (bt == 4L) "fixed" else if (bt == 5L) "none" else "gabriel")
  formats <- .cli_first(opt, c("formats"), NULL)
  if (is.null(formats)) {
    formats <- c(if (isTRUE(opt$OutPdf)) "pdf", if (isTRUE(opt$OutPng)) "png")
    if (!length(formats)) formats <- c("pdf", "svg")
  } else formats <- strsplit(formats, ",", fixed = TRUE)[[1L]]
  block_args <- list()
  if (!is.null(opt$BlockCut)) {
    z <- as.numeric(strsplit(opt$BlockCut, ":", fixed = TRUE)[[1L]])
    if (length(z) == 2L) {
      block_args$strong_cut <- z[1L]
      block_args$strong_fraction <- z[2L]
    }
  }
  if (!is.null(opt$FixBlock)) {
    block_args$fixed <- utils::read.table(opt$FixBlock, header = FALSE,
                                          col.names = c("chr", "start", "end"),
                                          stringsAsFactors = FALSE)
  }
  special <- NULL
  special_file <- .cli_first(opt, c("SpeSNPName", "special"))
  if (!is.null(special_file)) {
    special <- utils::read.table(special_file, header = FALSE,
                                 col.names = c("chr", "pos", "label"),
                                 stringsAsFactors = FALSE)
  }
  gwas <- .cli_first(opt, c("gwas", "InGWAS"))
  if (!is.null(gwas) && isTRUE(opt$NoLogP)) gwas <- read_gwas(gwas, region, value_is_logp = TRUE)
  result <- ldblockr(
    input = input, output_prefix = out, region = region, format = format,
    samples = .cli_first(opt, c("samples", "SubPop")),
    min_maf = as.numeric(.cli_first(opt, c("maf", "MAF"), 0.05)),
    max_missing = as.numeric(.cli_first(opt, c("missing", "Miss"), 0.25)),
    max_het = as.numeric(.cli_first(opt, c("het", "Het"), 1)),
    measure = measure,
    max_distance = as.numeric(.cli_first(opt, c("max-distance", "NoShowLDist"), Inf)),
    block_method = block_method,
    tag_threshold = as.numeric(.cli_first(opt, c("tag-r2", "TagSNPCut"), 0.8)),
    gwas = gwas, gff = .cli_first(opt, c("gff", "InGFF")), special = special,
    output_formats = formats,
    plot_args = list(cutline = as.numeric(.cli_first(opt, c("cutline", "Cutline"), Inf)),
                     lead = .cli_first(opt, c("lead", "TopSite"))),
    block_args = block_args
  )
  cat("LDblockR completed.\n")
  cat(paste(unname(result$files), collapse = "\n"), "\n")
  invisible(result)
}
