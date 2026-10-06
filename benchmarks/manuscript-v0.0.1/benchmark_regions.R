# One independent benchmark process. The wrapper records process peak RSS.
# Rscript benchmark_regions.R number_of_markers replicate output_directory [number_of_samples]
suppressPackageStartupMessages(library(LDblockR))
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3L || length(args) > 4L) stop("Expected marker count, replicate, output directory, and optional sample count")
m <- as.integer(args[1]); replicate <- as.integer(args[2]); output <- args[3]
n <- if(length(args) == 4L) as.integer(args[4]) else 281L
if(anyNA(c(m, replicate, n)) || m < 2L || n < 2L) stop("Marker and sample counts must be integers >= 2")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
set.seed(20260902 + m + n + replicate)
freq<-runif(m,.05,.5)
g<-vapply(freq,function(p)rbinom(n,2,p),numeric(n))
g[matrix(runif(n*m)<.01,n,m)]<-NA_real_
x<-as_ld_data(g,data.frame(chr="1",pos=seq_len(m)*1000,id=paste0("v",seq_len(m))))
gc()
calculation<-system.time(ld<-ld_compute(x,"r2",r2_method="dosage",min_n=20))[["elapsed"]]
p<-plot_ld(ld,"r2",show_maf=FALSE,position_scale="index",render="auto",draw=FALSE)
file<-file.path(output,paste0("benchmark_",m,"_",n,"_",replicate,".svg"))
rendering<-system.time(save_ld_plot(p,file,width=7.1,height=4.2))[["elapsed"]]
tab<-data.frame(n_samples=n,n_markers=m,replicate=replicate,ld_seconds=calculation,
  render_seconds=rendering,svg_bytes=file.info(file)$size,
  ld_object_bytes=as.numeric(object.size(ld)),
  rendering=if(m<=250)"vector" else "raster",R_version=R.version.string,package_version=as.character(packageVersion("LDblockR")))
write.table(tab,file.path(output,paste0("benchmark_",m,"_",n,"_",replicate,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
