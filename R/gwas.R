# Fast mixed-model association scan
#
# The implementation follows the same statistical model used by the
# P3D/MLM family of GAPIT analyses, but keeps the implementation independent:
# the null variance ratio is estimated once, the kinship matrix is
# eigendecomposed once, and SNPs are tested in BLAS-friendly chunks.

.gwas_as_path <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || !nzchar(x) || !file.exists(x)) {
    .stopf("%s must be an existing file path.", what)
  }
  normalizePath(x, mustWork = FALSE)
}

#' Read a TASSEL/rTASSEL phenotype table
#'
#' @param file A TASSEL phenotype file. The optional `<Phenotype>` and type
#'   rows are skipped automatically.
#' @return A data frame with the `Taxa` identifier column and numeric traits.
#' @export
read_tassel_phenotype <- function(file) {
  file <- .gwas_as_path(file, "file")
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) .stopf("Phenotype file is empty: %s", file)
  header <- which(grepl("(^|\\t)Taxa(\\t|$)", lines))[1L]
  if (!is.finite(header)) .stopf("Could not locate a Taxa header in %s.", file)
  tab <- utils::read.table(text = paste(lines[header:length(lines)], collapse = "\n"),
                           header = TRUE, sep = "\t", quote = "",
                           comment.char = "", fill = TRUE,
                           na.strings = c("", "NA", "NaN", "."),
                           check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) .stopf("No phenotype rows were found in %s.", file)
  id_col <- .match_column(tab, c("Taxa", "taxa", "sample", "id"), TRUE, "Taxa")
  names(tab)[match(id_col, names(tab))] <- "Taxa"
  for (j in seq_along(tab)) {
    if (j == match("Taxa", names(tab))) next
    z <- suppressWarnings(as.numeric(as.character(tab[[j]])))
    raw <- trimws(as.character(tab[[j]]))
    convertible <- any(is.finite(z)) && all(is.na(z) | raw == "" | is.finite(z))
    if (convertible) tab[[j]] <- z
  }
  # Repeated Taxa rows are valid in the TASSEL example (one row per
  # location).  Keep duplicates so `replicate = "expand"` can model them.
  tab$Taxa <- as.character(tab$Taxa)
  attr(tab, "source") <- file
  class(tab) <- c("tassel_phenotype", class(tab))
  tab
}

#' Read a TASSEL text kinship matrix
#'
#' @param file A TASSEL kinship text file.
#' @param samples Optional sample IDs to retain, in the desired order.
#' @return A symmetric numeric matrix with taxa as row and column names.
#' @export
read_tassel_kinship <- function(file, samples = NULL) {
  file <- .gwas_as_path(file, "file")
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, header = FALSE, sep = "\t", quote = "",
                           comment.char = "", fill = TRUE, check.names = FALSE,
                           na.strings = c("", "NA", "NaN", "."),
                           stringsAsFactors = FALSE)
  if (nrow(tab) < 2L || ncol(tab) < 3L) .stopf("Malformed TASSEL kinship file: %s", file)
  # TASSEL writes a leading dimension row (for example, `277` followed by
  # empty cells).  Data rows always start with a taxon name.
  first_nonempty <- trimws(as.character(tab[1L, ]))
  # `nzchar(NA)` is TRUE in base R, so remove missing padding before testing
  # whether the row contains only the leading dimension token.
  first_nonempty <- first_nonempty[!is.na(first_nonempty) & nzchar(first_nonempty) &
                                     !toupper(first_nonempty) %in% c("NA", "NAN", ".")]
  if (length(first_nonempty) == 1L && grepl("^[0-9]+$", first_nonempty)) tab <- tab[-1L, , drop = FALSE]
  ids <- as.character(tab[[1L]])
  vals <- tab[, -1L, drop = FALSE]
  k <- matrix(suppressWarnings(as.numeric(unlist(vals, use.names = FALSE))),
              nrow = nrow(vals), ncol = ncol(vals), byrow = FALSE)
  if (nrow(k) != ncol(k)) .stopf("Kinship matrix is not square in %s.", file)
  if (anyNA(k)) .stopf("Kinship matrix contains non-numeric values in %s.", file)
  ids <- make.unique(ids)
  rownames(k) <- ids
  colnames(k) <- ids
  k <- (k + t(k)) / 2
  if (!is.null(samples)) {
    samples <- .read_sample_ids(samples)
    idx <- match(samples, ids)
    if (anyNA(idx)) .stopf("Samples not found in kinship matrix: %s", paste(samples[is.na(idx)], collapse = ", "))
    k <- k[idx, idx, drop = FALSE]
  }
  attr(k, "source") <- file
  k
}

#' Return paths to built-in demonstration data
#'
#' @param name Data bundle name: `"tassel"` or `"regional"`.
#' @return A named character vector of installed example-file paths.
#' @export
example_data <- function(name = "tassel") {
  name <- match.arg(tolower(name), c("tassel", "regional"))
  if (name == "regional") {
    root <- system.file("extdata", package = "LDblockR")
    if (!nzchar(root)) .stopf("The built-in regional data bundle is not installed.")
    file_names <- c("example.vcf.gz", "example.hmp.txt", "example_gwas.tsv",
                    "regional_gwas.tsv", "example.gff3", "example_fixed_blocks.tsv",
                    "example_special.tsv", "subpopulation_A.txt", "subpopulation_B.txt",
                    "sample_groups.tsv", "example_matrix.tsv", "example_map.tsv",
                    "example_regions.tsv")
    keys <- c("vcf", "hapmap", "example_gwas", "regional_gwas", "gff3",
              "fixed_blocks", "special", "subpopulation_A", "subpopulation_B",
              "sample_groups", "matrix", "map", "regions")
    files <- file.path(root, file_names)
    plink_prefix <- file.path(root, "plink", "example_plink")
    files <- c(files, plink = plink_prefix)
    keys <- c(keys, "plink")
    required <- c(file.path(root, file_names), paste0(plink_prefix, c(".bed", ".bim", ".fam")))
    if (any(!file.exists(required))) {
      missing <- required[!file.exists(required)]
      .stopf("The built-in regional data bundle is incomplete; missing: %s. Reinstall the full LDblockR source package.",
             paste(basename(missing), collapse = ", "))
    }
    return(stats::setNames(files, keys))
  }
  root <- system.file("extdata", name, package = "LDblockR")
  if (!nzchar(root)) .stopf("The built-in %s data bundle is not installed.", name)
  files <- file.path(root, c("mdp_genotype.hmp.txt", "mdp_phenotype.txt",
                             "mdp_kinship.txt", "mdp_population_structure.txt",
                             "mdp_traits.txt", "SOURCE.md", "PROVENANCE.json"))
  if (any(!file.exists(files))) {
    .stopf("The built-in %s data bundle is incomplete; missing: %s. Reinstall the full LDblockR source package.",
           name, paste(basename(files[!file.exists(files)]), collapse = ", "))
  }
  stats::setNames(files, c("mdp_genotype", "mdp_phenotype", "mdp_kinship",
                           "mdp_population_structure", "mdp_traits", "SOURCE",
                           "PROVENANCE"))
}

.gwas_read_genotype <- function(genotype, map = NULL, sample_ids = NULL) {
  if (inherits(genotype, "ld_data")) {
    return(list(genotypes = genotype$genotypes, map = genotype$variants,
                samples = genotype$samples, source = genotype$source))
  }
  if (is.character(genotype) && length(genotype) == 1L) {
    file <- .gwas_as_path(genotype, "genotype")
    low <- tolower(file)
    if (grepl("(?:hmp|hapmap)(?:\\.txt)?(?:\\.gz)?$", low)) {
      x <- read_hapmap(file, min_maf = 0, max_missing = 1, max_het = 1,
                       quiet = TRUE)
      return(list(genotypes = x$genotypes, map = x$variants,
                  samples = x$samples, source = file))
    }
    if (grepl("\\.rds$", low)) {
      x <- readRDS(file)
      return(.gwas_read_genotype(x, map = map, sample_ids = sample_ids))
    }
    .stopf("Unsupported genotype file for gwas_mlm: %s. Use HapMap or an ld_data object.", file)
  }
  g <- as.matrix(genotype)
  storage.mode(g) <- "numeric"
  if (!nrow(g) || !ncol(g)) .stopf("genotype must be a non-empty samples-by-markers matrix.")
  bad <- !is.na(g) & !(g %in% c(0, 1, 2))
  if (any(bad)) .stopf("genotype must be coded 0, 1, 2, or NA.")
  ids <- sample_ids
  if (is.null(ids)) ids <- rownames(g)
  if (is.null(ids)) ids <- paste0("sample", seq_len(nrow(g)))
  if (length(ids) != nrow(g)) .stopf("sample_ids does not match genotype rows.")
  ids <- make.unique(as.character(ids))
  rownames(g) <- ids
  variants <- .normalize_map(map, ncol(g), colnames(g))
  colnames(g) <- variants$id
  list(genotypes = g, map = variants, samples = ids, source = NULL)
}

.gwas_numeric_matrix <- function(x, ids, label = "covariates") {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  if (!is.null(rownames(x)) && all(ids %in% rownames(x))) x <- x[ids, , drop = FALSE]
  if (nrow(x) != length(ids)) .stopf("%s has %d rows but phenotype has %d rows.", label, nrow(x), length(ids))
  for (j in seq_along(x)) {
    if (is.character(x[[j]])) {
      z <- suppressWarnings(as.numeric(x[[j]]))
      if (all(is.na(x[[j]]) | is.finite(z))) x[[j]] <- z
    }
  }
  x
}

.gwas_aggregate_phenotype <- function(pheno, id, y, covar, mode) {
  if (mode == "expand") return(list(id = id, y = y, covar = covar))
  groups <- split(seq_along(id), id, drop = TRUE)
  ids <- names(groups)
  yy <- vapply(groups, function(ii) mean(y[ii], na.rm = TRUE), numeric(1L))
  if (!is.null(covar)) {
    cv0 <- as.data.frame(covar, stringsAsFactors = FALSE)
    cv <- cv0[match(ids, id), , drop = FALSE]
    for (j in seq_along(cv0)) {
      if (is.numeric(cv0[[j]])) {
        cv[[j]] <- vapply(groups, function(ii) mean(cv0[[j]][ii], na.rm = TRUE), numeric(1L))
      } else {
        cv[[j]] <- vapply(groups, function(ii) as.character(cv0[[j]][ii[1L]]), character(1L))
      }
    }
  } else cv <- NULL
  list(id = ids, y = yy, covar = cv)
}

.gwas_relationship <- function(g) {
  means <- colMeans(g, na.rm = TRUE)
  for (j in seq_len(ncol(g))) {
    miss <- is.na(g[, j])
    if (any(miss)) g[miss, j] <- means[j]
  }
  z <- sweep(g, 2L, means, "-")
  s <- sqrt(colSums(z^2) / pmax(1, nrow(z) - 1L))
  keep <- is.finite(s) & s > 0
  if (!any(keep)) return(diag(nrow(g)))
  z <- sweep(z[, keep, drop = FALSE], 2L, s[keep], "/")
  k <- tcrossprod(z) / ncol(z)
  k <- (k + t(k)) / 2
  diag(k) <- pmax(diag(k), 1e-8)
  k
}

.gwas_independent_columns <- function(x) {
  q <- qr(x, tol = 1e-10)
  if (q$rank == ncol(x)) return(x)
  .warnf("Dropped %d linearly dependent fixed-effect columns.", ncol(x) - q$rank)
  x[, q$pivot[seq_len(q$rank)], drop = FALSE]
}

.gwas_profile_delta <- function(y, x, kinship, max_iter = 60L) {
  n <- length(y)
  p <- ncol(x)
  ev <- eigen((kinship + t(kinship)) / 2, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  if (!any(lam > 0)) lam <- rep(1, length(lam))
  yt <- drop(crossprod(ev$vectors, y))
  xt <- crossprod(ev$vectors, x)
  objective <- function(log_delta) {
    delta <- exp(log_delta)
    w <- 1 / sqrt(pmax(lam + delta, 1e-12))
    xw <- xt * w
    yw <- yt * w
    fit <- tryCatch(qr(xw, tol = 1e-10), error = function(e) NULL)
    if (is.null(fit) || fit$rank < p) return(Inf)
    beta <- qr.coef(fit, yw)
    res <- yw - drop(xw %*% beta)
    sse <- sum(res^2)
    if (!is.finite(sse) || sse <= 0) return(Inf)
    ld_x <- tryCatch(as.numeric(determinant(crossprod(xw), logarithm = TRUE)$modulus),
                     error = function(e) NA_real_)
    if (!is.finite(ld_x)) return(Inf)
    0.5 * ((n - p) * log(sse / max(1, n - p)) +
             sum(log(pmax(lam + delta, 1e-12))) + ld_x)
  }
  opt <- optimize(objective, interval = c(-8, 8), maximum = FALSE,
                  tol = 1e-5)
  delta <- exp(opt$minimum)
  list(delta = delta, eigen = ev, objective = opt$objective,
       log_delta = opt$minimum, iterations = max_iter)
}

.gwas_fit_null <- function(y, x, kinship, max_iter = 60L) {
  prof <- .gwas_profile_delta(y, x, kinship, max_iter = max_iter)
  ev <- prof$eigen
  lam <- pmax(ev$values, 0)
  delta <- prof$delta
  yt <- drop(crossprod(ev$vectors, y))
  xt <- crossprod(ev$vectors, x)
  w <- 1 / sqrt(pmax(lam + delta, 1e-12))
  xw <- xt * w
  yw <- yt * w
  fit <- qr(xw, tol = 1e-10)
  beta <- qr.coef(fit, yw)
  q <- qr.Q(fit, complete = FALSE)[, seq_len(fit$rank), drop = FALSE]
  residual <- yw - drop(xw %*% beta)
  sse <- sum(residual^2)
  df <- max(1, length(y) - fit$rank)
  list(delta = delta, eigen = ev, lambda = lam, transform = w,
       x_eigen = xt, x_white = xw, y_white = yw, q = q,
       beta = beta, residual = residual, sse = sse,
       df = df, rank = fit$rank, pve = 1 / (1 + delta),
       loglik_profile = -prof$objective)
}

.gwas_resolve_covariates <- function(phenotype, covariates, id, n_pc, g_taxa, taxa_ids) {
  if (is.null(covariates) && is.data.frame(phenotype)) {
    preferred <- intersect(c("location", "Q1", "Q2", "Q3"), names(phenotype))
    if (length(preferred)) covariates <- preferred
  }
  if (is.character(covariates) && is.data.frame(phenotype)) {
    miss <- setdiff(covariates, names(phenotype))
    if (length(miss)) .stopf("Covariate columns not found: %s", paste(miss, collapse = ", "))
    cv <- phenotype[match(id, phenotype$Taxa), covariates, drop = FALSE]
  } else cv <- .gwas_numeric_matrix(covariates, id, "covariates")
  if (is.null(cv)) cv <- data.frame(row.names = seq_along(id))
  if (n_pc > 0L) {
    z <- g_taxa
    means <- colMeans(z, na.rm = TRUE)
    for (j in seq_len(ncol(z))) {
      miss <- is.na(z[, j])
      if (any(miss)) z[miss, j] <- means[j]
    }
    z <- sweep(z, 2L, colMeans(z), "-")
    z <- sweep(z, 2L, apply(z, 2L, function(a) max(sd(a), 1e-12)), "/")
    sv <- svd(z, nu = min(nrow(z), n_pc), nv = 0L)
    pcs <- sv$u[, seq_len(min(n_pc, ncol(sv$u))), drop = FALSE]
    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
    pcs <- pcs[match(id, taxa_ids), , drop = FALSE]
    cv <- cbind(cv, as.data.frame(pcs, stringsAsFactors = FALSE))
  }
  cv
}

#' Fast GAPIT-principle mixed linear model GWAS
#'
#' Fits the model `y = X beta + Zu + e`, with `u ~ N(0, K sigma_g^2)`.
#' The null variance ratio is estimated once (P3D), followed by a vectorized
#' score/GLS scan in marker chunks. It is the fixed-variance single-marker MLM
#' test used by the P3D family while avoiding a new REML optimization for every
#' marker. `PVE` is the marker-specific partial phenotypic variance
#' explained; `model_PVE` is the null-model variance component estimate.
#'
#' @param genotype An `ld_data` object, a 0/1/2 matrix, a HapMap file, or an
#'   RDS file containing one of those objects.
#' @param phenotype A phenotype data frame, a TASSEL phenotype file, or a
#'   numeric vector with names matching genotype samples.
#' @param trait Trait column name. The first numeric phenotype column is used
#'   when omitted.
#' @param kinship Optional numeric kinship matrix or TASSEL kinship file. If
#'   omitted, a standardized marker relationship is computed once.
#' @param covariates Fixed-effect columns, a matrix/data frame, or `NULL` to
#'   use `location`, `Q1`, `Q2`, and `Q3` when present.
#' @param n_pc Number of genotype PCs to append to the fixed effects.
#' @param map Optional marker map when `genotype` is a matrix.
#' @param sample_ids Optional sample IDs when matrix row names are absent.
#' @param replicate How repeated phenotype rows are handled: `"expand"`
#'   keeps location replicates, while `"mean"` averages them by taxon.
#' @param min_maf Minimum minor allele frequency.
#' @param max_missing Maximum marker missing-call fraction.
#' @param chunk_size Number of markers processed per BLAS chunk.
#' @param max_iter Reserved compatibility argument for variance estimation.
#' @param verbose Print model and progress messages.
#' @return A data frame of association statistics with class
#'   `gwas_mlm_result`; model details are stored in `attr(x, "model")`.
#' @export
gwas_mlm <- function(genotype, phenotype, trait = NULL, kinship = NULL,
                     covariates = NULL, n_pc = 0L, map = NULL,
                     sample_ids = NULL, replicate = c("expand", "mean"),
                     min_maf = 0.05, max_missing = 0.2,
                     chunk_size = 256L, max_iter = 60L, verbose = TRUE) {
  replicate <- match.arg(replicate)
  chunk_size <- max(1L, as.integer(chunk_size))
  n_pc <- max(0L, as.integer(n_pc))
  g <- .gwas_read_genotype(genotype, map = map, sample_ids = sample_ids)
  gmat <- g$genotypes
  marker_map <- g$map
  taxa <- g$samples

  if (is.character(phenotype) && length(phenotype) == 1L) phenotype <- read_tassel_phenotype(phenotype)
  if (is.numeric(phenotype) && is.null(dim(phenotype))) {
    y <- as.numeric(phenotype)
    id <- names(phenotype)
    if (is.null(id)) id <- taxa[seq_along(y)]
    ph <- data.frame(Taxa = id, trait_value = y, stringsAsFactors = FALSE)
    trait <- "trait_value"
  } else {
    ph <- as.data.frame(phenotype, stringsAsFactors = FALSE)
    if (!nrow(ph)) .stopf("phenotype is empty.")
    id_col <- .match_column(ph, c("Taxa", "taxa", "sample", "id"), TRUE, "Taxa")
    ph$Taxa <- as.character(ph[[id_col]])
    if (is.null(trait)) {
      candidates <- setdiff(names(ph), c(id_col, "Taxa"))
      numeric_candidates <- candidates[vapply(ph[candidates], function(z) {
        zz <- suppressWarnings(as.numeric(as.character(z)))
        any(is.finite(zz)) && !all(is.na(zz))
      }, logical(1L))]
      if (!length(numeric_candidates)) .stopf("No numeric trait column was found in phenotype.")
      trait <- numeric_candidates[1L]
    }
    if (!trait %in% names(ph)) .stopf("Trait column '%s' was not found in phenotype.", trait)
    y <- suppressWarnings(as.numeric(as.character(ph[[trait]])))
    id <- ph$Taxa
  }
  if (length(y) != nrow(ph)) .stopf("Phenotype trait length is inconsistent.")
  keep <- is.finite(y) & !is.na(id) & nzchar(id)
  if (!any(keep)) .stopf("No finite phenotype observations remain.")
  ph <- ph[keep, , drop = FALSE]
  y <- y[keep]
  id <- as.character(ph$Taxa)
  taxa_idx <- match(id, taxa)
  keep <- !is.na(taxa_idx)
  if (!any(keep)) .stopf("No phenotype Taxa IDs overlap genotype samples.")
  y <- y[keep]
  id <- id[keep]
  taxa_idx <- taxa_idx[keep]
  ph <- ph[keep, , drop = FALSE]

  if (is.character(covariates) && is.data.frame(ph)) {
    miss <- setdiff(covariates, names(ph))
    if (length(miss)) .stopf("Covariate columns not found: %s", paste(miss, collapse = ", "))
    cv <- ph[, covariates, drop = FALSE]
  } else if (is.null(covariates)) {
    preferred <- intersect(c("location", "Q1", "Q2", "Q3"), names(ph))
    cv <- if (length(preferred)) ph[, preferred, drop = FALSE] else NULL
  } else cv <- covariates
  agg <- .gwas_aggregate_phenotype(ph, id, y, cv, replicate)
  id <- agg$id
  y <- agg$y
  cv <- agg$covar
  taxa_idx <- match(id, taxa)
  if (anyNA(taxa_idx)) .stopf("Phenotype aggregation produced unknown genotype IDs.")
  if (is.null(cv)) cv <- data.frame(row.names = seq_along(id))
  if (nrow(cv) != length(y)) .stopf("Covariate rows do not match phenotype observations.")

  # Build/align the taxa-level kinship before expanding repeated environments.
  if (is.null(kinship)) {
    k_taxa <- .gwas_relationship(gmat)
    rownames(k_taxa) <- colnames(k_taxa) <- taxa
  } else if (is.character(kinship) && length(kinship) == 1L) {
    k_taxa <- read_tassel_kinship(kinship)
    shared <- intersect(taxa, rownames(k_taxa))
    if (length(shared) < 3L) .stopf("Kinship and genotype have fewer than three shared taxa.")
    g_keep <- match(shared, taxa)
    gmat <- gmat[g_keep, , drop = FALSE]
    taxa <- taxa[g_keep]
    marker_map <- marker_map
    taxa_idx <- match(id, taxa)
    keep <- !is.na(taxa_idx)
    y <- y[keep]; id <- id[keep]; taxa_idx <- taxa_idx[keep]
    cv <- cv[keep, , drop = FALSE]
    k_taxa <- k_taxa[shared, shared, drop = FALSE]
  } else {
    k_taxa <- as.matrix(kinship)
    storage.mode(k_taxa) <- "numeric"
    if (is.null(colnames(k_taxa)) && !is.null(rownames(k_taxa))) {
      colnames(k_taxa) <- rownames(k_taxa)
    }
    if (is.null(rownames(k_taxa))) {
      if (nrow(k_taxa) != length(taxa)) .stopf("Unnamed kinship must have one row per genotype sample.")
      rownames(k_taxa) <- colnames(k_taxa) <- taxa
    }
    shared <- intersect(taxa, rownames(k_taxa))
    if (length(shared) < 3L) .stopf("Kinship and genotype have fewer than three shared taxa.")
    g_keep <- match(shared, taxa)
    gmat <- gmat[g_keep, , drop = FALSE]
    taxa <- taxa[g_keep]
    taxa_idx <- match(id, taxa)
    keep <- !is.na(taxa_idx)
    y <- y[keep]; id <- id[keep]; taxa_idx <- taxa_idx[keep]
    cv <- cv[keep, , drop = FALSE]
    k_taxa <- k_taxa[shared, shared, drop = FALSE]
  }
  if (length(y) < 10L) .stopf("At least ten overlapping phenotype observations are required.")
  k_obs <- k_taxa[taxa_idx, taxa_idx, drop = FALSE]
  rownames(k_obs) <- colnames(k_obs) <- make.unique(id)

  cv <- .gwas_resolve_covariates(ph, covariates, id, n_pc, gmat, taxa)
  # For aggregated phenotypes, the helper above uses the original `ph` rows;
  # use the already aggregated covariate table when it is available.
  if (!is.null(agg$covar)) cv <- agg$covar
  if (!nrow(cv)) cv <- data.frame(row.names = seq_along(y))
  x <- if (ncol(cv)) stats::model.matrix(~ ., data = cv) else matrix(1, nrow = length(y), ncol = 1L,
                                                                      dimnames = list(NULL, "(Intercept)"))
  x <- .gwas_independent_columns(x)
  if (nrow(x) != length(y)) .stopf("Fixed-effect design does not match phenotype observations.")

  model <- .gwas_fit_null(y, x, k_obs, max_iter = max_iter)
  if (verbose) {
    message(sprintf("MLM null model: n=%d, markers=%d, fixed effects=%d, delta=%.4g, model PVE=%.3f",
                    length(y), ncol(gmat), ncol(x), model$delta, model$pve))
  }

  obs_g <- gmat[taxa_idx, , drop = FALSE]
  called <- colSums(!is.na(gmat))
  af <- colMeans(gmat, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  missing_rate <- 1 - called / nrow(gmat)
  keep_markers <- is.finite(maf) & maf >= min_maf & missing_rate <= max_missing & called >= 3L
  if (!any(keep_markers)) .stopf("No markers pass min_maf=%g and max_missing=%g.", min_maf, max_missing)
  marker_idx <- which(keep_markers)
  result <- vector("list", ceiling(length(marker_idx) / chunk_size))
  n_done <- 0L
  for (chunk_start in seq(1L, length(marker_idx), by = chunk_size)) {
    ii <- marker_idx[chunk_start:min(length(marker_idx), chunk_start + chunk_size - 1L)]
    gc <- obs_g[, ii, drop = FALSE]
    # Mean imputation is performed only for the tested marker matrix.  The
    # reported MAF and missingness always come from the observed calls.
    for (j in seq_len(ncol(gc))) {
      miss <- is.na(gc[, j])
      if (any(miss)) gc[miss, j] <- mean(gc[, j], na.rm = TRUE)
    }
    ge <- crossprod(model$eigen$vectors, gc)
    gw <- sweep(ge, 1L, model$transform, "*")
    proj <- crossprod(model$q, gw)
    gr <- gw - model$q %*% proj
    denom <- colSums(gr^2)
    num <- drop(crossprod(model$residual, gr))
    valid <- is.finite(denom) & denom > 1e-12 & is.finite(num)
    beta <- rep(NA_real_, length(denom))
    beta[valid] <- num[valid] / denom[valid]
    sse_full <- rep(NA_real_, length(denom))
    sse_full[valid] <- pmax(0, model$sse - num[valid]^2 / denom[valid])
    df_full <- max(1, model$df - 1L)
    se <- rep(NA_real_, length(denom))
    se[valid] <- sqrt((sse_full[valid] / df_full) / denom[valid])
    stat <- beta / se
    p <- 2 * stats::pt(-abs(stat), df = df_full)
    p[!is.finite(p)] <- NA_real_
    pve <- rep(NA_real_, length(denom))
    pve[valid] <- pmax(0, pmin(1, (model$sse - sse_full[valid]) / model$sse))
    allele <- ifelse(!is.na(marker_map$ref[ii]) & !is.na(marker_map$alt[ii]),
                     paste0(marker_map$ref[ii], "/", marker_map$alt[ii]), "")
    result[[ceiling(chunk_start / chunk_size)]] <- data.frame(
      chr = as.character(marker_map$chr[ii]), pos = as.numeric(marker_map$pos[ii]),
      id = as.character(marker_map$id[ii]), allele = allele,
      maf = maf[ii], n_called = called[ii], missing_rate = missing_rate[ii],
      beta = beta, se = se, statistic = stat, p = p,
      logp = .safe_log10(p), PVE = pve, model_PVE = rep(model$pve, length(ii)),
      stringsAsFactors = FALSE
    )
    n_done <- n_done + length(ii)
    if (verbose && (n_done == length(marker_idx) || n_done %% (chunk_size * 4L) == 0L)) {
      message(sprintf("  tested %d/%d markers", n_done, length(marker_idx)))
    }
  }
  out <- do.call(rbind, result)
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$pos, out$id), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "model") <- list(trait = trait, n = length(y), n_taxa = length(unique(id)),
                              n_markers_tested = nrow(out), fixed_effects = colnames(x),
                              delta = model$delta, model_PVE = model$pve,
                              sigma_g2 = model$sse / model$df,
                              sigma_e2 = model$delta * model$sse / model$df,
                              replicate = replicate, kinship_source = attr(k_taxa, "source"),
                              genotype_source = g$source, min_maf = min_maf,
                              max_missing = max_missing)
  attr(out, "call") <- match.call()
  class(out) <- c("gwas_mlm_result", class(out))
  out
}

#' Select a regional LD window from threshold-exceeding GWAS markers
#'
#' The selected locus is restricted to one local cluster on one chromosome
#' (the cluster with the strongest peak, then the most hits) so that the
#' returned interval can be passed directly to `read_hapmap()`,
#' `read_vcf_region()`, or `read_plink()`. If no marker reaches the threshold,
#' the most significant marker is used as a reproducible fallback and
#' `threshold_hit` is set to `FALSE` in the returned table.
#'
#' @param gwas A table or file accepted by [read_gwas()].
#' @param cutline Threshold on `-log10(P)`. If `NULL`, use the Bonferroni
#'   threshold `-log10(0.05/nrow(gwas))`.
#' @param flank Number of bases added on each side of the selected hits.
#' @param min_width Minimum interval width in bases.
#' @param max_width Optional maximum interval width; it must still contain all
#'   threshold-exceeding hits on the selected chromosome.
#' @param chromosome Optional chromosome to force.
#' @return A list of class `gwas_region` with `region`, `chr`, `start`,
#'   `end`, `lead_id`, `hit_ids`, `cutline`, and the selected GWAS table.
#' @export
gwas_ld_region <- function(gwas, cutline = NULL, flank = 250000,
                           min_width = 0, max_width = Inf,
                           chromosome = NULL) {
  tab <- read_gwas(gwas)
  if (!nrow(tab)) .stopf("GWAS table has no finite rows.")
  if (is.null(cutline) || length(cutline) == 0L || is.na(cutline[1L])) {
    cutline <- -log10(0.05 / nrow(tab))
  }
  cutline <- as.numeric(cutline[1L])
  flank <- as.numeric(flank[1L])
  min_width <- as.numeric(min_width[1L])
  max_width <- as.numeric(max_width[1L])
  if (!is.finite(cutline) || cutline < 0) .stopf("cutline must be a finite non-negative number.")
  if (!is.finite(flank) || flank < 0) .stopf("flank must be a non-negative number.")
  if (!is.finite(min_width) || min_width < 0) .stopf("min_width must be a non-negative number.")
  if (is.na(max_width) || max_width <= 0) .stopf("max_width must be positive or Inf.")
  if (max_width < min_width) .stopf("max_width cannot be smaller than min_width.")
  hit <- is.finite(tab$logp) & tab$logp >= cutline
  available_chr <- unique(as.character(tab$chr))
  if (!is.null(chromosome)) {
    chromosome <- as.character(chromosome[1L])
    if (!chromosome %in% available_chr) .stopf("Chromosome '%s' is absent from the GWAS table.", chromosome)
    chr <- chromosome
  } else if (any(hit)) {
    by_chr <- table(as.character(tab$chr[hit]))
    candidates <- names(by_chr)[by_chr == max(by_chr)]
    if (length(candidates) == 1L) {
      chr <- candidates
    } else {
      score <- vapply(candidates, function(z) {
        max(tab$logp[hit & as.character(tab$chr) == z], na.rm = TRUE)
      }, numeric(1L))
      chr <- candidates[which.max(score)]
    }
  } else {
    chr <- as.character(tab$chr[which.max(tab$logp)])
  }
  chr_idx <- which(as.character(tab$chr) == chr)
  chr_hit_all <- chr_idx[hit[chr_idx]]
  # Separate distant significant loci on the same chromosome.  The LD panel
  # should describe one local locus rather than silently spanning unrelated
  # peaks; the default gap is twice the requested flank.
  if (length(chr_hit_all) > 1L) {
    ord <- order(tab$pos[chr_hit_all])
    sorted_hit <- chr_hit_all[ord]
    grp <- cumsum(c(TRUE, diff(tab$pos[sorted_hit]) > max(1, 2 * flank)))
    clusters <- split(sorted_hit, grp)
    score <- vapply(clusters, function(ii) max(tab$logp[ii], na.rm = TRUE), numeric(1L))
    count <- vapply(clusters, length, integer(1L))
    best <- order(-score, -count)[1L]
    chr_hit <- clusters[[best]]
  } else {
    chr_hit <- chr_hit_all
  }
  threshold_hit <- length(chr_hit) > 0L
  anchor_idx <- if (threshold_hit) chr_hit[which.max(tab$logp[chr_hit])] else {
    chr_idx[which.max(tab$logp[chr_idx])]
  }
  range_idx <- if (threshold_hit) chr_hit else anchor_idx
  hit_start <- min(tab$pos[range_idx], na.rm = TRUE)
  hit_end <- max(tab$pos[range_idx], na.rm = TRUE)
  hit_span <- max(0, hit_end - hit_start)
  requested_width <- max(hit_span + 2 * flank, min_width)
  if (is.finite(max_width) && max_width < hit_span) {
    .stopf("max_width (%.0f) is smaller than the span of threshold hits (%.0f).",
           max_width, hit_span)
  }
  width <- if (is.finite(max_width)) min(requested_width, max_width) else requested_width
  # Center on the hit span so every selected threshold marker remains inside
  # the returned interval, including when max_width caps the flank.
  center <- (hit_start + hit_end) / 2
  start <- floor(center - width / 2)
  end <- ceiling(center + width / 2)
  if (start < 0) {
    end <- end - start
    start <- 0
  }
  if (is.finite(max_width) && end - start > max_width) end <- start + max_width
  selected <- tab[as.character(tab$chr) == chr & tab$pos >= start & tab$pos <= end,
                  , drop = FALSE]
  selected$threshold_hit <- selected$logp >= cutline
  selected <- selected[order(selected$pos, selected$id), , drop = FALSE]
  rownames(selected) <- NULL
  lead_id <- as.character(tab$id[anchor_idx])
  out <- list(
    region = paste0(chr, ":", format(start, scientific = FALSE, trim = TRUE), "-",
                    format(end, scientific = FALSE, trim = TRUE)),
    chr = chr, start = start, end = end, cutline = cutline,
    flank = flank, lead_id = lead_id,
    hit_ids = as.character(tab$id[chr_hit]),
    threshold_hit = threshold_hit, selected = selected
  )
  class(out) <- c("gwas_region", "list")
  out
}

#' @export
print.gwas_region <- function(x, ...) {
  cat("<gwas_region>", x$region, "| cutline:", format(x$cutline, digits = 4), "\n")
  cat("  lead:", x$lead_id, "| threshold hits:", length(x$hit_ids), "\n")
  invisible(x)
}

#' @export
print.gwas_mlm_result <- function(x, ...) {
  model <- attr(x, "model")
  cat("<gwas_mlm_result>", nrow(x), "markers tested\n")
  if (!is.null(model)) {
    cat("  trait:", model$trait, "| observations:", model$n,
        "| model PVE:", format(model$model_PVE, digits = 3), "\n")
  }
  if (nrow(x)) {
    # Subsetting a gwas_mlm_result retains its S3 class.  Strip that class
    # before printing the preview, otherwise print.gwas_mlm_result() calls
    # itself recursively through print.data.frame and eventually overflows
    # the C stack on ordinary head()/subset() workflows.
    cols <- intersect(c("id", "chr", "pos", "p", "PVE"), names(x))
    top <- x[order(x$p, na.last = TRUE), cols, drop = FALSE]
    class(top) <- setdiff(class(top), "gwas_mlm_result")
    print(utils::head(top, 5L), row.names = FALSE)
  }
  invisible(x)
}
