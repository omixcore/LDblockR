set.seed(507)

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_samples <- 60L
n_variants <- 42L
samples <- sprintf("Line%03d", seq_len(n_samples))
positions <- 1000000 + cumsum(sample(900:2100, n_variants, replace = TRUE))
ids <- paste0("snp", seq_len(n_variants))
ref <- rep(c("A", "C", "G", "T"), length.out = n_variants)
alt <- c(A = "G", C = "T", G = "A", T = "C")[ref]

hap <- matrix(0L, nrow = 2L * n_samples, ncol = n_variants)
blocks <- list(1:14, 15:28, 29:42)
block_freq <- c(0.28, 0.42, 0.20)
for (b in seq_along(blocks)) {
  latent <- rbinom(2L * n_samples, 1L, block_freq[b])
  for (j in blocks[[b]]) {
    mutation <- rbinom(2L * n_samples, 1L, if (j %% 7L == 0L) 0.12 else 0.025)
    hap[, j] <- abs(latent - mutation)
  }
}
geno <- hap[seq(1, 2L * n_samples, by = 2L), ] + hap[seq(2, 2L * n_samples, by = 2L), ]
missing <- matrix(runif(length(geno)) < 0.008, nrow = n_samples)
geno[missing] <- NA_integer_

vcf <- c(
  "##fileformat=VCFv4.2",
  "##source=LDblockR-example",
  "##contig=<ID=chr1,length=2000000>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t")
)
for (j in seq_len(n_variants)) {
  gt <- character(n_samples)
  for (i in seq_len(n_samples)) {
    if (is.na(geno[i, j])) gt[i] <- ".|."
    else gt[i] <- paste0(hap[2L * i - 1L, j], "|", hap[2L * i, j])
  }
  vcf <- c(vcf, paste(c("chr1", positions[j], ids[j], ref[j], alt[j], 60, "PASS", ".", "GT", gt), collapse = "\t"))
}
con <- gzfile(file.path(out_dir, "example.vcf.gz"), "wt")
writeLines(vcf, con)
close(con)

r2 <- suppressWarnings(cor(geno, use = "pairwise.complete.obs")^2)
lead <- 21L
logp <- 1.2 + 8.5 * r2[lead, ] + runif(n_variants, 0, 0.65)
logp[lead] <- 12.5
gwas <- data.frame(chr = "chr1", pos = positions, p = 10^(-logp), id = ids)
write.table(gwas, file.path(out_dir, "example_gwas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Keep the documented regional-workflow filename in sync with the generated
# association table.  The package also resolves this bare name after install.
file.copy(file.path(out_dir, "example_gwas.tsv"),
          file.path(out_dir, "regional_gwas.tsv"), overwrite = TRUE)

gff <- c(
  "##gff-version 3",
  paste("chr1", "LDblockR", "gene", positions[2] - 1500, positions[13] + 1200, ".", "+", ".", "ID=ZmGene01;Name=ZmGene01", sep = "\t"),
  paste("chr1", "LDblockR", "mRNA", positions[2] - 1500, positions[13] + 1200, ".", "+", ".", "ID=ZmGene01_T01;Parent=ZmGene01;Name=ZmGene01_T01", sep = "\t"),
  paste("chr1", "LDblockR", "five_prime_UTR", positions[2] - 1500, positions[3], ".", "+", ".", "ID=utr1;Parent=ZmGene01_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[3], positions[6], ".", "+", "0", "ID=cds1;Parent=ZmGene01_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[9], positions[12], ".", "+", "0", "ID=cds2;Parent=ZmGene01_T01", sep = "\t"),
  paste("chr1", "LDblockR", "three_prime_UTR", positions[12], positions[13] + 1200, ".", "+", ".", "ID=utr2;Parent=ZmGene01_T01", sep = "\t"),
  paste("chr1", "LDblockR", "gene", positions[17] - 800, positions[28] + 600, ".", "-", ".", "ID=ZmGene02;Name=ZmGene02", sep = "\t"),
  paste("chr1", "LDblockR", "mRNA", positions[17] - 800, positions[28] + 600, ".", "-", ".", "ID=ZmGene02_T01;Parent=ZmGene02;Name=ZmGene02_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[18], positions[21], ".", "-", "0", "ID=cds3;Parent=ZmGene02_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[24], positions[27], ".", "-", "0", "ID=cds4;Parent=ZmGene02_T01", sep = "\t"),
  paste("chr1", "LDblockR", "gene", positions[31] - 1000, positions[40] + 1000, ".", "+", ".", "ID=ZmGene03;Name=ZmGene03", sep = "\t"),
  paste("chr1", "LDblockR", "mRNA", positions[31] - 1000, positions[40] + 1000, ".", "+", ".", "ID=ZmGene03_T01;Parent=ZmGene03;Name=ZmGene03_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[32], positions[35], ".", "+", "0", "ID=cds5;Parent=ZmGene03_T01", sep = "\t"),
  paste("chr1", "LDblockR", "CDS", positions[37], positions[40], ".", "+", "0", "ID=cds6;Parent=ZmGene03_T01", sep = "\t")
)
writeLines(gff, file.path(out_dir, "example.gff3"))

hapmap <- data.frame(
  `rs#` = ids, alleles = paste0(ref, "/", alt), chrom = "chr1", pos = positions,
  strand = "+", `assembly#` = "v4", center = "LDblockR", protLSID = "NA",
  assayLSID = "NA", panelLSID = "NA", QCcode = "QC+", check.names = FALSE
)
base_calls <- matrix("N", nrow = n_variants, ncol = n_samples)
for (j in seq_len(n_variants)) {
  for (i in seq_len(n_samples)) {
    if (is.na(geno[i, j])) base_calls[j, i] <- "N"
    else if (geno[i, j] == 0) base_calls[j, i] <- ref[j]
    else if (geno[i, j] == 2) base_calls[j, i] <- alt[j]
    else base_calls[j, i] <- paste0(ref[j], alt[j])
  }
}
colnames(base_calls) <- samples
hapmap <- cbind(hapmap, as.data.frame(base_calls, check.names = FALSE))
write.table(hapmap, file.path(out_dir, "example.hmp.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

matrix_out <- data.frame(sample = samples, geno, check.names = FALSE)
names(matrix_out)[-1L] <- ids
write.table(matrix_out, file.path(out_dir, "example_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(chr = "chr1", pos = positions, id = ids, ref = ref, alt = alt),
            file.path(out_dir, "example_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(samples[seq_len(30L)], file.path(out_dir, "subpopulation_A.txt"))
write.table(data.frame(chr = "chr1", start = c(positions[1], positions[15], positions[29]),
                       end = c(positions[14], positions[28], positions[42])),
            file.path(out_dir, "example_fixed_blocks.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(data.frame(chr = "chr1", pos = positions[c(7, 21, 35)],
                       label = c("Candidate-A", "Lead-SNP", "Candidate-B")),
            file.path(out_dir, "example_special.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(data.frame(chr = "chr1", start = c(positions[1], positions[15]),
                       end = c(positions[20], positions[42])),
            file.path(out_dir, "example_regions.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
