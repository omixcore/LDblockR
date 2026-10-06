.extract_gt <- function(sample_fields, gt_index) {
  vapply(sample_fields, function(z) {
    p <- strsplit(z, ":", fixed = TRUE)[[1L]]
    if (length(p) >= gt_index) p[gt_index] else "."
  }, character(1L), USE.NAMES = FALSE)
}

.decode_vcf_gt <- function(gt) {
  missing <- is.na(gt) | gt %in% c(".", "./.", ".|.")
  has_sep <- grepl("[/|]", gt)
  a1s <- sub("[/|].*$", "", gt)
  a2s <- ifelse(has_sep, sub("^.*[/|]", "", gt), a1s)
  a1 <- suppressWarnings(as.integer(a1s))
  a2 <- suppressWarnings(as.integer(a2s))
  invalid <- missing | is.na(a1) | is.na(a2) | a1 < 0L | a1 > 1L | a2 < 0L | a2 > 1L
  dosage <- a1 + a2
  dosage[invalid] <- NA_real_
  phased <- grepl("\\|", gt)
  unphased_het <- !invalid & dosage == 1 & !phased
  h1 <- as.numeric(a1)
  h2 <- as.numeric(a2)
  h1[invalid | unphased_het] <- NA_real_
  h2[invalid | unphased_het] <- NA_real_
  list(dosage = dosage, h1 = h1, h2 = h2,
       phase_complete = !any(unphased_het, na.rm = TRUE))
}

.read_vcf_base <- function(file, region = NULL, samples = NULL,
                           pass_only = FALSE, biallelic_snps_only = TRUE,
                           chunk_size = 5000L) {
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  header <- NULL
  repeat {
    z <- readLines(con, n = 1L, warn = FALSE)
    if (!length(z)) break
    if (startsWith(z, "#CHROM")) {
      header <- strsplit(sub("^#", "", z), "\t", fixed = TRUE)[[1L]]
      break
    }
  }
  if (is.null(header) || length(header) < 10L) .stopf("VCF header with sample columns was not found in %s.", file)
  all_samples <- header[10L:length(header)]
  requested <- .read_sample_ids(samples)
  if (is.null(requested)) {
    sample_index <- seq_along(all_samples)
    selected_samples <- all_samples
  } else {
    sample_index <- match(requested, all_samples)
    if (anyNA(sample_index)) .stopf("Samples not found in VCF: %s", paste(requested[is.na(sample_index)], collapse = ", "))
    selected_samples <- requested
  }
  geno_list <- list()
  hap_list <- list()
  map_list <- list()
  k <- 0L
  done <- FALSE
  while (!done) {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    fields <- strsplit(lines, "\t", fixed = TRUE)
    for (f in fields) {
      if (length(f) < 9L) next
      chr <- f[1L]
      pos <- suppressWarnings(as.numeric(f[2L]))
      if (!is.finite(pos)) next
      if (!is.null(region)) {
        if (!.same_chr(chr, region$chr[1L]) || pos < region$start[1L]) next
        if (.same_chr(chr, region$chr[1L]) && pos > region$end[1L]) {
          done <- TRUE
          break
        }
      }
      if (pass_only && !(f[7L] %in% c("PASS", "."))) next
      if (biallelic_snps_only && (grepl(",", f[5L], fixed = TRUE) || nchar(f[4L]) != 1L || nchar(f[5L]) != 1L)) next
      fmt <- strsplit(f[9L], ":", fixed = TRUE)[[1L]]
      gt_index <- match("GT", fmt)
      if (is.na(gt_index)) next
      sf <- f[9L + sample_index]
      gt <- .extract_gt(sf, gt_index)
      dec <- .decode_vcf_gt(gt)
      k <- k + 1L
      geno_list[[k]] <- dec$dosage
      h <- as.vector(rbind(dec$h1, dec$h2))
      hap_list[[k]] <- h
      id <- f[3L]
      if (is.na(id) || id == "." || !nzchar(id)) id <- paste0(chr, ":", format(pos, scientific = FALSE, trim = TRUE))
      map_list[[k]] <- data.frame(
        chr = chr, pos = pos, id = id, ref = f[4L], alt = f[5L],
        qual = suppressWarnings(as.numeric(f[6L])), filter = f[7L],
        phase_complete = dec$phase_complete, stringsAsFactors = FALSE
      )
    }
  }
  if (!k) .stopf("No usable variants were found in %s%s.", file,
                  if (is.null(region)) "" else paste0(" for ", region$chr, ":", region$start, "-", region$end))
  geno <- do.call(cbind, geno_list)
  hap <- do.call(cbind, hap_list)
  map <- do.call(rbind, map_list)
  as_ld_data(geno, map, selected_samples, haplotypes = hap, source = normalizePath(file, mustWork = FALSE))
}

.bcftools_region_extract <- function(file, region) {
  exe <- Sys.which("bcftools")
  if (!nzchar(exe) || is.null(region)) return(NULL)
  tmp <- tempfile(fileext = ".vcf")
  err <- tempfile(fileext = ".log")
  reg <- paste0(region$chr[1L], ":", region$start[1L], "-", region$end[1L])
  status <- suppressWarnings(system2(exe, c("view", "-r", shQuote(reg), shQuote(file)), stdout = tmp, stderr = err))
  if (!identical(status, 0L)) {
    unlink(tmp)
    unlink(err)
    return(NULL)
  }
  unlink(err)
  tmp
}

#' Read a VCF region
#'
#' Reads plain or gzip-compressed VCF, supports sample subsets, and uses
#' `bcftools` automatically when available for indexed regional extraction.
#'
#' @param file VCF or VCF.GZ file.
#' @param region Optional `chr:start-end` string.
#' @param samples Optional character vector or one-column sample-list file.
#' @param min_maf Minimum minor allele frequency.
#' @param max_missing Maximum missing fraction.
#' @param max_het Maximum heterozygous-call fraction.
#' @param min_hwe_p Minimum HWE P value; zero disables this filter.
#' @param pass_only Keep only PASS or unfiltered records.
#' @param biallelic_snps_only Exclude indels and multiallelic records.
#' @param backend `auto`, `base`, or `bcftools`.
#' @param max_variants Maximum variants to retain.
#' @param on_excess Error or evenly thin when `max_variants` is exceeded.
#' @param quiet Suppress filtering messages.
#' @return An `ld_data` object.
#' @export
read_vcf_region <- function(file, region = NULL, samples = NULL,
                            min_maf = 0.05, max_missing = 0.25,
                            max_het = 1, min_hwe_p = 0,
                            pass_only = FALSE, biallelic_snps_only = TRUE,
                            backend = c("auto", "base", "bcftools"),
                            max_variants = Inf, on_excess = c("error", "thin"),
                            quiet = FALSE) {
  backend <- match.arg(backend)
  on_excess <- match.arg(on_excess)
  reg <- parse_region(region)
  extracted <- NULL
  use_bcftools <- backend == "bcftools" || (backend == "auto" && !is.null(reg) && nzchar(Sys.which("bcftools")))
  if (use_bcftools) {
    extracted <- .bcftools_region_extract(file, reg)
    if (is.null(extracted) && backend == "bcftools") .stopf("bcftools regional extraction failed; check that the VCF is indexed and the region exists.")
    if (is.null(extracted) && backend == "auto") .warnf("bcftools extraction failed; falling back to the built-in streaming VCF reader.")
  }
  parse_file <- if (is.null(extracted)) file else extracted
  on.exit(if (!is.null(extracted)) unlink(extracted), add = TRUE)
  x <- .read_vcf_base(parse_file, if (is.null(extracted)) reg else NULL, samples,
                      pass_only, biallelic_snps_only)
  x$source <- normalizePath(file, mustWork = FALSE)
  x <- filter_variants(x, min_maf = min_maf, max_missing = max_missing,
                       max_het = max_het, min_hwe_p = min_hwe_p, quiet = quiet)
  if (x$n_variants > max_variants) {
    if (on_excess == "error") {
      .stopf("%d variants remain, exceeding max_variants=%d. Narrow the region, raise max_variants, or use on_excess='thin'.",
             x$n_variants, as.integer(max_variants))
    }
    keep <- unique(round(seq(1, x$n_variants, length.out = as.integer(max_variants))))
    x <- subset_ld_data(x, variants = keep)
    .warnf("Variants were evenly thinned to %d for regional analysis.", x$n_variants)
  }
  x
}

.decode_hapmap_value <- function(z, ref, alt) {
  z <- toupper(gsub("[/|[:space:]]", "", as.character(z)))
  if (!nzchar(z) || z %in% c("N", "NN", "-", "--", ".")) return(NA_real_)
  iupac <- c(R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
  if (nchar(z) == 1L && z %in% names(iupac)) z <- iupac[[z]]
  if (nchar(z) == 1L) z <- paste0(z, z)
  if (nchar(z) != 2L) return(NA_real_)
  a <- strsplit(z, "", fixed = TRUE)[[1L]]
  if (any(!a %in% c(ref, alt))) return(NA_real_)
  sum(a == alt)
}

#' Read HapMap genotype data
#'
#' @param file HapMap text file, optionally compressed.
#' @param region Optional genomic region.
#' @param samples Optional sample vector or list file.
#' @param min_maf,max_missing,max_het,min_hwe_p Variant filters.
#' @param quiet Suppress filtering messages.
#' @return An `ld_data` object.
#' @export
read_hapmap <- function(file, region = NULL, samples = NULL,
                        min_maf = 0.05, max_missing = 0.25,
                        max_het = 1, min_hwe_p = 0, quiet = FALSE) {
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, header = TRUE, sep = "\t", quote = "",
                           comment.char = "", check.names = FALSE,
                           stringsAsFactors = FALSE)
  if (ncol(tab) < 12L) .stopf("HapMap input must contain the 11 standard metadata columns plus samples.")
  id_col <- .match_column(tab, c("rs#", "id", "snp"), TRUE, "marker ID")
  allele_col <- .match_column(tab, c("alleles", "allele"), TRUE, "alleles")
  chr_col <- .match_column(tab, c("chrom", "chr", "chromosome"), TRUE, "chromosome")
  pos_col <- .match_column(tab, c("pos", "position", "bp"), TRUE, "position")
  reg <- parse_region(region)
  keep_rows <- .in_region(tab[[chr_col]], as.numeric(tab[[pos_col]]), reg)
  tab <- tab[keep_rows, , drop = FALSE]
  if (!nrow(tab)) .stopf("No HapMap variants occur in the requested region.")
  sample_cols <- names(tab)[12L:ncol(tab)]
  requested <- .read_sample_ids(samples)
  if (!is.null(requested)) {
    miss <- setdiff(requested, sample_cols)
    if (length(miss)) .stopf("Samples not found in HapMap file: %s", paste(miss, collapse = ", "))
    sample_cols <- requested
  }
  allele_parts <- strsplit(as.character(tab[[allele_col]]), "[/|]", perl = TRUE)
  ref <- vapply(allele_parts, function(z) if (length(z)) toupper(z[1L]) else NA_character_, character(1L))
  alt <- vapply(allele_parts, function(z) if (length(z) >= 2L) toupper(z[2L]) else NA_character_, character(1L))
  valid <- !is.na(ref) & !is.na(alt) & nchar(ref) == 1L & nchar(alt) == 1L
  tab <- tab[valid, , drop = FALSE]
  ref <- ref[valid]
  alt <- alt[valid]
  if (!nrow(tab)) .stopf("No biallelic HapMap variants were found.")
  geno <- matrix(NA_real_, nrow = length(sample_cols), ncol = nrow(tab),
                 dimnames = list(sample_cols, as.character(tab[[id_col]])))
  for (j in seq_len(nrow(tab))) {
    geno[, j] <- vapply(tab[j, sample_cols, drop = TRUE], .decode_hapmap_value,
                        numeric(1L), ref = ref[j], alt = alt[j])
  }
  map <- data.frame(chr = as.character(tab[[chr_col]]), pos = as.numeric(tab[[pos_col]]),
                    id = as.character(tab[[id_col]]), ref = ref, alt = alt,
                    stringsAsFactors = FALSE)
  x <- as_ld_data(geno, map, sample_cols, source = normalizePath(file, mustWork = FALSE))
  filter_variants(x, min_maf = min_maf, max_missing = max_missing,
                  max_het = max_het, min_hwe_p = min_hwe_p, quiet = quiet)
}

.read_plink_bed <- function(prefix) {
  bed <- paste0(prefix, ".bed")
  bim <- paste0(prefix, ".bim")
  fam <- paste0(prefix, ".fam")
  if (!all(file.exists(c(bed, bim, fam)))) .stopf("PLINK BED input requires %s.bed, .bim, and .fam.", prefix)
  b <- utils::read.table(bim, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  f <- utils::read.table(fam, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  if (ncol(b) < 6L || ncol(f) < 2L) .stopf("Malformed PLINK BIM or FAM file.")
  n <- nrow(f)
  m <- nrow(b)
  con <- file(bed, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 3L)
  if (length(magic) != 3L || !identical(as.integer(magic), c(108L, 27L, 1L))) {
    .stopf("Unsupported PLINK BED header. SNP-major BED format is required.")
  }
  bytes_per_variant <- ceiling(n / 4)
  raw <- readBin(con, what = "raw", n = bytes_per_variant * m)
  if (length(raw) != bytes_per_variant * m) .stopf("PLINK BED file is truncated.")
  geno <- matrix(NA_real_, nrow = n, ncol = m)
  code_to_dosage <- c(0, NA, 1, 2)
  raw_int <- as.integer(raw)
  for (j in seq_len(m)) {
    offset <- (j - 1L) * bytes_per_variant
    for (i in seq_len(n)) {
      byte <- raw_int[offset + ((i - 1L) %/% 4L) + 1L]
      code <- bitwAnd(bitwShiftR(byte, 2L * ((i - 1L) %% 4L)), 3L)
      geno[i, j] <- code_to_dosage[code + 1L]
    }
  }
  ids <- make.unique(as.character(f[[2L]]))
  rownames(geno) <- ids
  colnames(geno) <- as.character(b[[2L]])
  map <- data.frame(chr = as.character(b[[1L]]), id = as.character(b[[2L]]),
                    cm = as.numeric(b[[3L]]), pos = as.numeric(b[[4L]]),
                    ref = as.character(b[[5L]]), alt = as.character(b[[6L]]),
                    stringsAsFactors = FALSE)
  list(genotypes = geno, map = map, samples = ids)
}

.read_plink_ped <- function(prefix) {
  ped_file <- paste0(prefix, ".ped")
  map_file <- paste0(prefix, ".map")
  if (!all(file.exists(c(ped_file, map_file)))) .stopf("PLINK PED input requires %s.ped and .map.", prefix)
  mp <- utils::read.table(map_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  ped <- utils::read.table(ped_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  m <- nrow(mp)
  if (ncol(mp) < 4L || ncol(ped) != 6L + 2L * m) .stopf("Malformed PLINK PED/MAP input.")
  geno <- matrix(NA_real_, nrow = nrow(ped), ncol = m)
  ref <- alt <- rep(NA_character_, m)
  for (j in seq_len(m)) {
    a1 <- as.character(ped[[6L + 2L * j - 1L]])
    a2 <- as.character(ped[[6L + 2L * j]])
    called <- c(a1, a2)
    called <- called[!called %in% c("0", "N", ".", "-")]
    alleles <- unique(called)
    if (length(alleles) < 2L) alleles <- c(alleles, NA_character_)
    ref[j] <- alleles[1L]
    alt[j] <- alleles[2L]
    miss <- a1 %in% c("0", "N", ".", "-") | a2 %in% c("0", "N", ".", "-") | is.na(alt[j])
    geno[, j] <- as.numeric(a1 == alt[j]) + as.numeric(a2 == alt[j])
    geno[miss, j] <- NA_real_
  }
  samples <- make.unique(as.character(ped[[2L]]))
  rownames(geno) <- samples
  colnames(geno) <- as.character(mp[[2L]])
  map <- data.frame(chr = as.character(mp[[1L]]), id = as.character(mp[[2L]]),
                    cm = as.numeric(mp[[3L]]), pos = as.numeric(mp[[4L]]),
                    ref = ref, alt = alt, stringsAsFactors = FALSE)
  list(genotypes = geno, map = map, samples = samples)
}

#' Read PLINK genotype files
#'
#' @param prefix Prefix for BED/BIM/FAM or PED/MAP files.
#' @param region Optional region.
#' @param samples Optional sample vector or list file.
#' @param min_maf,max_missing,max_het,min_hwe_p Variant filters.
#' @param quiet Suppress filtering messages.
#' @return An `ld_data` object.
#' @export
read_plink <- function(prefix, region = NULL, samples = NULL,
                       min_maf = 0.05, max_missing = 0.25,
                       max_het = 1, min_hwe_p = 0, quiet = FALSE) {
  raw <- if (file.exists(paste0(prefix, ".bed"))) .read_plink_bed(prefix) else .read_plink_ped(prefix)
  x <- as_ld_data(raw$genotypes, raw$map, raw$samples,
                  source = normalizePath(prefix, mustWork = FALSE))
  reg <- parse_region(region)
  if (!is.null(reg)) {
    keep <- .in_region(x$variants$chr, x$variants$pos, reg)
    if (!any(keep)) .stopf("No PLINK variants occur in the requested region.")
    x <- subset_ld_data(x, variants = keep)
  }
  requested <- .read_sample_ids(samples)
  if (!is.null(requested)) x <- subset_ld_data(x, samples = requested)
  filter_variants(x, min_maf = min_maf, max_missing = max_missing,
                  max_het = max_het, min_hwe_p = min_hwe_p, quiet = quiet)
}

#' Read genotype data with automatic format detection
#'
#' @param x File path, PLINK prefix, matrix, or existing `ld_data` object.
#' @param format `auto`, `vcf`, `hapmap`, `plink`, or `matrix`.
#' @param ... Arguments forwarded to the format-specific reader. For matrix
#'   input, `map`, `sample_ids`, `haplotypes`, and `source` are accepted.
#' @return An `ld_data` object.
#' @export
read_genotypes <- function(x, format = c("auto", "vcf", "hapmap", "plink", "matrix"), ...) {
  format <- match.arg(format)
  if (inherits(x, "ld_data")) return(x)
  if (is.matrix(x) || is.data.frame(x)) {
    if (format == "auto") format <- "matrix"
    if (format != "matrix") .stopf("In-memory table input requires format='matrix'.")
    return(as_ld_data(x, ...))
  }
  if (!is.character(x) || length(x) != 1L) .stopf("x must be a file path, PLINK prefix, matrix, or ld_data object.")
  if (format == "auto") {
    low <- tolower(x)
    if (grepl("\\.vcf(?:\\.gz|\\.bgz)?$", low)) format <- "vcf"
    else if (grepl("(?:hmp|hapmap)(?:\\.txt)?(?:\\.gz)?$", low)) format <- "hapmap"
    else if (file.exists(paste0(x, ".bed")) || file.exists(paste0(x, ".ped"))) format <- "plink"
    else format <- "matrix"
  }
  if (format == "vcf") return(read_vcf_region(x, ...))
  if (format == "hapmap") return(read_hapmap(x, ...))
  if (format == "plink") return(read_plink(x, ...))
  args <- list(...)
  map <- args$map
  sample_ids <- args$sample_ids
  haplotypes <- args$haplotypes
  source <- .null_coalesce(args$source, normalizePath(x, mustWork = FALSE))
  tab <- utils::read.table(x, header = TRUE, row.names = 1L, check.names = FALSE,
                           comment.char = "", stringsAsFactors = FALSE)
  as_ld_data(as.matrix(tab), map = map, sample_ids = sample_ids,
             haplotypes = haplotypes, source = source)
}
