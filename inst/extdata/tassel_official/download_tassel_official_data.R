# Download full official TASSEL/rTASSEL example files
# Run from the root of LDblockR_demo_complete.
out_dir <- file.path("data", "tassel_official")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  maize_chr9_10thin40000.recode.vcf = "https://raw.githubusercontent.com/maize-genetics/rTASSEL/master/inst/extdata/maize_chr9_10thin40000.recode.vcf",
  mdp_genotype.hmp.txt = "https://raw.githubusercontent.com/maize-genetics/rTASSEL/master/inst/extdata/mdp_genotype.hmp.txt",
  mdp_phenotype.txt = "https://raw.githubusercontent.com/maize-genetics/rTASSEL/master/inst/extdata/mdp_phenotype.txt",
  mdp_population_structure.txt = "https://raw.githubusercontent.com/maize-genetics/rTASSEL/master/inst/extdata/mdp_population_structure.txt",
  wheat_gbs_test.vcf = "https://raw.githubusercontent.com/maize-genetics/rTASSEL/master/inst/extdata/wheat_gbs_test.vcf"
)

for (nm in names(files)) {
  dest <- file.path(out_dir, nm)
  message("Downloading: ", nm)
  tryCatch(
    download.file(files[[nm]], destfile = dest, mode = "wb", quiet = FALSE),
    error = function(e) message("  FAILED: ", conditionMessage(e))
  )
}
message("Finished. Files are in: ", normalizePath(out_dir, mustWork = FALSE))
