# LDblockR 0.0.1 manuscript reproduction. Only R and LDblockR are required.
# Rscript reproduce_manuscript.R output_directory
# Or source(system.file("scripts", "reproduce_manuscript.R", package="LDblockR"))
#    result <- reproduce_ldblockr("LDblockR_results")

example_files <- function(name = "regional") {
  name <- match.arg(name, c("regional", "tassel"))
  root <- system.file("extdata", package = "LDblockR")
  if (!nzchar(root)) stop("The installed LDblockR package has no extdata directory.")
  if (name == "regional") {
    out <- c(
      vcf = file.path(root, "example.vcf.gz"),
      hapmap = file.path(root, "example.hmp.txt"),
      matrix = file.path(root, "example_matrix.tsv"),
      map = file.path(root, "example_map.tsv"),
      plink = file.path(root, "plink", "example_plink"),
      gwas = file.path(root, "regional_gwas.tsv"),
      gff = file.path(root, "example.gff3")
    )
  } else {
    tassel_root <- file.path(root, "tassel")
    out <- c(
      genotype = file.path(tassel_root, "mdp_genotype.hmp.txt"),
      phenotype = file.path(tassel_root, "mdp_phenotype.txt"),
      kinship = file.path(tassel_root, "mdp_kinship.txt")
    )
  }
  required <- out[!grepl("plink$", names(out))]
  if ("plink" %in% names(out)) {
    required <- c(required, paste0(out[["plink"]], c(".bed", ".bim", ".fam",
                                                        ".ped", ".map")))
  }
  if (any(!file.exists(required))) {
    stop("Missing installed example files: ", paste(basename(required[!file.exists(required)]), collapse = ", "))
  }
  out
}

reproduce_ldblockr <- function(output_dir = "LDblockR_results", figures = TRUE) {
  if (!requireNamespace("LDblockR", quietly = TRUE)) stop("Install LDblockR first.")
  if (utils::packageVersion("LDblockR") < "0.0.1") stop("LDblockR >= 0.0.1 is required.")
  suppressPackageStartupMessages(library(LDblockR))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  for (d in c("results", "figures", "projects")) dir.create(file.path(output_dir,d), showWarnings=FALSE)
  write_tab <- function(x, name) utils::write.table(x,file.path(output_dir,"results",name),
    sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  reference <- function(name) system.file("extdata","reference",name,package="LDblockR")
  ref_pair <- function(ld, tab, metric) ld[[metric]][cbind(match(tab$variant_1,ld$data$variants$id),
                                                         match(tab$variant_2,ld$data$variants$id))]
  validation <- list()
  add_validation <- function(test, n, difference, tolerance) {
    ok <- length(difference)==n && all(is.finite(difference)) && max(abs(difference))<=tolerance
    validation[[length(validation)+1L]] <<- data.frame(test=test,n=n,
      max_absolute_difference=max(abs(difference)),mean_absolute_difference=mean(abs(difference)),
      tolerance=tolerance,status=if(ok) "PASS" else "FAIL")
    if (!ok) stop("Validation failed: ", test)
  }
  message("1/4 Importing and comparing six synthetic input encodings")
  f <- example_files()
  synthetic <- read_vcf_region(f["vcf"],min_maf=0,max_missing=1,backend="base",quiet=TRUE)
  ped <- tempfile("ldblockr_ped_")
  on.exit(unlink(paste0(ped,c(".ped",".map"))),add=TRUE)
  file.copy(paste0(f["plink"],c(".ped",".map")),paste0(ped,c(".ped",".map")))
  inputs <- list(
    VCF=read_vcf_region(system.file("extdata","example.vcf",package="LDblockR"),min_maf=0,max_missing=1,quiet=TRUE),
    VCF_GZ=synthetic,
    HapMap=read_hapmap(f["hapmap"],min_maf=0,max_missing=1,quiet=TRUE),
    PLINK_BED=read_plink(f["plink"],min_maf=0,max_missing=1,quiet=TRUE),
    PLINK_PED=read_plink(ped,min_maf=0,max_missing=1,quiet=TRUE),
    Matrix=read_genotypes(f["matrix"],format="matrix",map=utils::read.table(f["map"],header=TRUE)))
  format_audit <- do.call(rbind,lapply(names(inputs),function(name) {
    x<-inputs[[name]]
    stopifnot(identical(x$samples,synthetic$samples),identical(x$variants$id,synthetic$variants$id),
              identical(x$variants$pos,synthetic$variants$pos))
    g<-x$genotypes
    flip<-which(!is.na(x$variants$ref) & x$variants$ref==synthetic$variants$alt &
                x$variants$alt==synthetic$variants$ref)
    if(length(flip))g[,flip]<-2-g[,flip]
    mask_mismatch<-sum(is.na(g)!=is.na(synthetic$genotypes))
    called_mismatch<-sum(g!=synthetic$genotypes,na.rm=TRUE)
    stopifnot(mask_mismatch==0,called_mismatch==0)
    data.frame(format=name,samples=x$n_samples,variants=x$n_variants,
      called_genotypes=sum(!is.na(g)),missing_genotypes=sum(is.na(g)),
      genotype_mismatches_vs_VCF=called_mismatch,missingness_mismatches_vs_VCF=mask_mismatch)
  }))
  write_tab(format_audit,"R_format_validation.tsv")
  phased <- ld_compute(synthetic,"both",r2_method="haplotype",min_n=5)
  dosage <- ld_compute(synthetic,"r2",r2_method="dosage",min_n=5)
  ref <- utils::read.delim(reference("synthetic_reference_pairs.tsv"),check.names=FALSE)
  ref$R_r2_phased <- ref_pair(phased,ref,"r2")
  ref$R_r2_dosage <- ref_pair(dosage,ref,"r2")
  ref$R_Dprime_phased <- ref_pair(phased,ref,"dprime")
  ref$R_n <- ref_pair(phased,ref,"n")
  stopifnot(nrow(ref)==861,all(ref$R_n==ref$n_diploid))
  add_validation("Synthetic phased r2 vs independent formula",861,ref$R_r2_phased-ref$r2_phased,1e-10)
  add_validation("Synthetic dosage r2 vs independent formula",861,ref$R_r2_dosage-ref$r2_dosage,1e-10)
  add_validation("Synthetic phased Dprime vs independent formula",861,ref$R_Dprime_phased-ref$abs_Dprime_phased,1e-10)
  message("2/4 Checking the archived, externally executed LDBlockShow outputs")
  triangle <- function(file,n) {
    con<-gzfile(file,"rt");on.exit(close(con))
    lines<-readLines(con,warn=FALSE)
    stopifnot(length(lines)==n,startsWith(lines[1],"#"))
    m<-matrix(NA_real_,n,n)
    for(i in seq_len(n-1L)) {
      z<-scan(text=lines[i+1L],quiet=TRUE)
      stopifnot(length(z)==n-i)
      m[i,(i+1L):n]<-z;m[(i+1L):n,i]<-z
    }
    m
  }
  extdir <- reference("ldblockshow")
  sites<-utils::read.table(gzfile(file.path(extdir,"synthetic.site.gz")),header=FALSE)
  stopifnot(identical(as.character(sites[[1]]),synthetic$variants$chr),
            all(sites[[2]]==synthetic$variants$pos))
  er2<-triangle(file.path(extdir,"synthetic.TriangleB.gz"),42)
  edp<-triangle(file.path(extdir,"synthetic.TriangleV.gz"),42)
  ij<-cbind(match(ref$variant_1,synthetic$variants$id),match(ref$variant_2,synthetic$variants$id))
  ref$LDBlockShow_r2_3dp<-er2[ij];ref$LDBlockShow_Dprime_3dp<-edp[ij]
  ref$r2_external_absolute_difference<-abs(ref$R_r2_phased-ref$LDBlockShow_r2_3dp)
  ref$Dprime_external_absolute_difference<-abs(ref$R_Dprime_phased-ref$LDBlockShow_Dprime_3dp)
  add_validation("Phased r2 vs LDBlockShow three-decimal output",861,ref$r2_external_absolute_difference,.00051)
  add_validation("Phased Dprime vs LDBlockShow three-decimal output",861,ref$Dprime_external_absolute_difference,.00051)
  write_tab(ref,"R_synthetic_pair_validation.tsv")
  message("3/4 Reproducing QC and the objectively selected TASSEL maize region")
  tf <- example_files("tassel")
  consistent <- read_hapmap(tf["genotype"],min_maf=0,max_missing=1,allele_conflict="drop_variant",quiet=TRUE)
  complete <- filter_variants(consistent,min_maf=0,max_missing=.2,quiet=TRUE)
  maize <- filter_variants(complete,min_maf=.05,max_missing=.2,max_het=1,min_hwe_p=0,quiet=TRUE)
  qc <- data.frame(step=c("Input markers","Declared-allele consistent","Missing fraction <= 0.20","MAF >= 0.05"),
                   retained_markers=c(nrow(consistent$input_qc),consistent$n_variants,complete$n_variants,maize$n_variants),
                   removed_markers=c(0,nrow(consistent$input_qc)-consistent$n_variants,
                                     consistent$n_variants-complete$n_variants,complete$n_variants-maize$n_variants))
  stopifnot(maize$n_samples==281,identical(qc$retained_markers,c(3093L,3045L,3045L,2526L)))
  write_tab(qc,"R_maize_qc_flow.tsv");write_tab(consistent$input_qc,"R_maize_import_decisions.tsv")
  write_tab(tail(maize$qc_filters,1)[[1]],"R_maize_final_filter_decisions.tsv")
  bychr<-data.frame(chromosome=as.character(1:10),raw_markers=as.integer(table(factor(consistent$input_qc$chr,levels=as.character(1:10)))),
    allele_consistent=as.integer(table(factor(consistent$variants$chr,levels=as.character(1:10)))),
    retained_markers=as.integer(table(factor(maize$variants$chr,levels=as.character(1:10)))))
  write_tab(bychr,"R_maize_chromosome_summary.tsv")
  v<-maize$variants
  counts<-vapply(seq_len(nrow(v)),function(i) sum(v$chr==v$chr[i] & v$pos>=v$pos[i] & v$pos<=v$pos[i]+999999),integer(1))
  first<-which.max(counts) # object is already ordered by chromosome, position and ID
  reg<-paste0(v$chr[first],":",format(v$pos[first],scientific=FALSE,trim=TRUE),"-",
              format(v$pos[first]+999999,scientific=FALSE,trim=TRUE))
  local<-subset_ld_data(maize,variants=which(v$chr==v$chr[first] & v$pos>=v$pos[first] & v$pos<=v$pos[first]+999999))
  maize_ld<-ld_compute(local,"r2",r2_method="dosage",min_n=20)
  mref<-utils::read.delim(reference("maize_regional_pairs.tsv"),check.names=FALSE)
  mref$R_r2_dosage<-ref_pair(maize_ld,mref,"r2");mref$R_n<-ref_pair(maize_ld,mref,"n")
  stopifnot(local$n_variants==16,all(mref$R_n==mref$n_diploid))
  add_validation("Maize regional dosage r2 vs independent formula",120,mref$R_r2_dosage-mref$r2_dosage,1e-10)
  write_tab(mref,"R_maize_pair_validation.tsv");write_tab(local$variants,"R_maize_regional_variants.tsv")
  summary<-data.frame(quantity=c("package_version","R_version","platform","maize_samples","maize_input_markers",
    "maize_raw_missing_calls","maize_conflicting_markers","maize_retained_markers","maize_region_AGPv1",
    "maize_regional_markers","maize_regional_pairs","maize_pairwise_n_min","maize_pairwise_n_max",
    "maize_median_r2","maize_mean_r2","maize_pairs_r2_at_least_0.8","maize_retained_median_MAF",
    "dosage_phased_mean_absolute_difference","dosage_phased_max_absolute_difference"),
    value=as.character(c(as.character(utils::packageVersion("LDblockR")),R.version.string,R.version$platform,281,3093,
      sum(consistent$input_qc$original_missing_calls),sum(consistent$input_qc$conflicting_calls>0),2526,reg,16,120,
      min(mref$R_n),max(mref$R_n),median(mref$R_r2_dosage),mean(mref$R_r2_dosage),sum(mref$R_r2_dosage>=.8),
      median(maize$variants$maf),mean(abs(ref$R_r2_dosage-ref$R_r2_phased)),max(abs(ref$R_r2_dosage-ref$R_r2_phased)))))
  write_tab(summary,"R_reproduction_summary.tsv")
  validation<-do.call(rbind,validation);write_tab(validation,"R_numerical_validation.tsv")
  blocks<-detect_ld_blocks(phased,"strong",metric="r2",strong_cut=.7,strong_fraction=.7)
  tags<-select_tag_snps(phased,threshold=.8,blocks=blocks)
  synthetic_plot<-plot_ld(phased,"r2",gwas=f["gwas"],genes=f["gff"],blocks=blocks,tags=tags,
    lead="snp21",cutline=5,title="Synthetic demonstration: 60 samples, 42 SNPs",palette="red",draw=FALSE)
  maize_plot<-plot_ld(maize_ld,"r2",position_scale="index",palette="red",title=paste("TASSEL maize tutorial |",reg,"(AGPv1)"),draw=FALSE)
  save_ld_project(phased,file.path(output_dir,"projects","synthetic_phased.rds"))
  save_ld_project(maize_ld,file.path(output_dir,"projects","maize_region.rds"))
  save_ld_project(synthetic_plot,file.path(output_dir,"projects","synthetic_plot.rds"))
  export_ld(phased,file.path(output_dir,"results","synthetic"),blocks,tags)
  export_ld(maize_ld,file.path(output_dir,"results","maize_region"))
  capture.output(utils::sessionInfo(),file=file.path(output_dir,"results","R_sessionInfo.txt"))
  result<-list(synthetic=synthetic,phased=phased,dosage=dosage,synthetic_pairs=ref,
    synthetic_plot=synthetic_plot,maize=maize,maize_ld=maize_ld,maize_plot=maize_plot,
    maize_region=reg,maize_pairs=mref,qc=qc,bychr=bychr,validation=validation,summary=summary,
    output_dir=output_dir)
  if(figures) {
    message("4/4 Drawing manuscript figures and standalone regional plots")
    env<-new.env(parent=globalenv())
    sys.source(system.file("scripts","publication_figures.R",package="LDblockR"),envir=env)
    env$make_publication_figures(result)
    for(ext in c("pdf","svg","png")) {
      save_ld_plot(synthetic_plot,file.path(output_dir,"figures",paste0("Synthetic_integrated_region.",ext)),width=7.1,height=5.3,dpi=600)
      save_ld_plot(maize_plot,file.path(output_dir,"figures",paste0("TASSEL_maize_region.",ext)),width=7.1,height=4.3,dpi=600)
    }
  }
  message("Completed: ",output_dir,"; all six numerical comparisons passed.")
  invisible(result)
}

if (sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  reproduce_ldblockr(if(length(args))args[1L] else "LDblockR_results")
}
