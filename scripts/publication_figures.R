# Publication figures constructed directly from executed LDblockR results.
# Edit this script to change typography, colors, or panel placement.

make_publication_figures <- function(result) {
  if (!requireNamespace("LDblockR", quietly = TRUE)) stop("Install LDblockR first.")
  stopifnot(is.list(result),all(result$validation$status=="PASS"))
  output <- file.path(result$output_dir,"figures")
  ink <- "#25313A"; muted <- "#58666D"; rule <- "#CDD5DA"
  blue <- "#2879A8"; teal <- "#23897C"; purple <- "#7961A5"; orange <- "#CA792F"
  page <- function() {
    graphics::par(mar=c(0,0,0,0),xaxs="i",yaxs="i",family="sans",fg=ink,col=ink,cex=1)
    graphics::plot.new();graphics::plot.window(c(0,1),c(0,1))
  }
  text <- function(x,y,label,cex=1,lineheight=1,...) {
    previous<-graphics::par("lheight");on.exit(graphics::par(lheight=previous))
    graphics::par(lheight=lineheight)
    graphics::text(x,y,label,cex=cex,col=ink,...)
  }
  label <- function(x,y,letter,title) {
    text(x,y,letter,cex=1.36,font=2,adj=c(0,1))
    text(x+.032,y-.002,title,cex=1.09,font=2,adj=c(0,1))
  }
  arrow <- function(x0,y0,x1,y1,col=muted) graphics::arrows(x0,y0,x1,y1,length=.065,lwd=.75,col=col)
  boxtext <- function(x0,x1,y0,y1,title,detail,color) {
    graphics::rect(x0,y0,x1,y1,col=grDevices::adjustcolor(color,alpha.f=.055),border=rule,lwd=.55)
    graphics::segments(x0,y0,x0,y1,col=color,lwd=2.2)
    text(x0+.011,y1-.013,title,cex=1,font=2,adj=c(0,1))
    text(x0+.011,y1-.042,detail,cex=.93,adj=c(0,1),lineheight=1.13)
  }
  figure1 <- function() {
    page()
    text(.02,.986,"LDblockR: reproducible regional LD analysis",cex=1.37,font=2,adj=c(0,1))
    text(.98,.981,"R package 0.0.1",cex=.96,adj=c(1,1))
    graphics::segments(.02,.945,.98,.945,col=rule,lwd=.7)
    label(.02,.927,"a","Input and quality control")
    boxtext(.025,.287,.797,.881,"VCF / VCF.GZ","Diploid calls; optional GT phase + PS",blue)
    boxtext(.025,.287,.704,.788,"HapMap / PLINK / matrix","IUPAC; BED or PED; 0/1/2/NA",teal)
    text(.026,.678,"MAF  |  missingness  |  heterozygosity",cex=.93,adj=c(0,1))
    text(.026,.653,"Optional HWE; recorded allele conflicts",cex=.93,adj=c(0,1))
    text(.026,.616,"Official TASSEL maize tutorial bundled",cex=.91,font=3,adj=c(0,1))
    arrow(.295,.762,.327,.762)
    label(.333,.927,"b","Traceable R objects")
    boxtext(.337,.650,.765,.881,"ld_data","genotypes; variants; samples\nhaplotypes; phase_sets; provenance",purple)
    arrow(.491,.760,.491,.741)
    boxtext(.337,.650,.646,.736,"ld_result","LD matrices; pairwise n; method flags\nparameters; analysis_history",blue)
    text(.341,.616,"Workflow project: data + results + plot",cex=.94,adj=c(0,1))
    arrow(.657,.762,.690,.762)
    label(.695,.927,"c","Analysis modules")
    text(.708,.872,"Pairwise LD",cex=1.03,font=2,adj=c(0,1))
    text(.708,.847,"Dosage or phased r²; absolute D′",cex=.94,adj=c(0,1))
    text(.708,.802,"Block definitions",cex=1.03,font=2,adj=c(0,1))
    text(.708,.777,"Strong-pair; solid-spine; four-gamete\nGabriel-inspired; supplied intervals",cex=.94,adj=c(0,1),lineheight=1.1)
    text(.708,.714,"Derived summaries",cex=1.03,font=2,adj=c(0,1))
    text(.708,.689,"Greedy tags; neighbors; LD decay\nGroup contrasts; haplotype frequencies",cex=.94,adj=c(0,1),lineheight=1.1)
    graphics::segments(.02,.579,.98,.579,col=rule,lwd=.7)
    label(.02,.560,"d","Integrated regional visualization")
    text(.058,.527,"Synthetic fixture: 60 samples, 42 SNPs; simulated association and genes",cex=.91,adj=c(0,1))
    label(.714,.560,"e","Numerical validation")
    text(.722,.505,"LDBlockShow 1.41",cex=1.03,font=2,adj=c(0,1))
    text(.722,.475,"861 / 861 pairs within 0.00051\nfor phased r² and absolute D′",cex=.96,adj=c(0,1),lineheight=1.15)
    v<-result$validation[4:5,]
    old<-graphics::par(no.readonly=TRUE)
    graphics::par(fig=c(.716,.985,.238,.430),mar=c(2.0,3.35,.8,.3),new=TRUE,cex=.88)
    values<-rbind(v$mean_absolute_difference,v$max_absolute_difference)*1e4
    graphics::barplot(values,beside=TRUE,col=c(blue,orange),border=NA,ylim=c(0,5.7),
      names.arg=c(expression(r^2),expression("|"*D*"'|")),las=1,cex.names=1,
      cex.axis=.86,ylab="",space=c(.2,.9))
    graphics::abline(h=5.1,lty=2,lwd=.85,col=muted)
    graphics::mtext(expression("Absolute difference ("*10^{-4}*")"),side=2,line=2.2,cex=.8)
    graphics::par(old)
    graphics::par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=TRUE);graphics::plot.new();graphics::plot.window(c(0,1),c(0,1))
    graphics::rect(.73,.221,.741,.231,col=blue,border=NA);text(.75,.226,"Mean",cex=.91,adj=0)
    graphics::rect(.83,.221,.841,.231,col=orange,border=NA);text(.85,.226,"Maximum",cex=.91,adj=0)
    text(.722,.195,"Dashed line: rounding tolerance.\nReference output has three decimals.",cex=.87,adj=c(0,1),lineheight=1.1)
    graphics::segments(.722,.135,.977,.135,col=rule,lwd=.6)
    text(.722,.117,"Reproducible outputs",cex=1.02,font=2,adj=c(0,1))
    text(.722,.088,"PDF / SVG / PNG / TIFF\nPair tables + saved RDS project\nInput hashes + parameters + R session",cex=.94,adj=c(0,1),lineheight=1.23)
    p<-result$synthetic_plot
    p$title<-NULL;p$viewport<-c(.017,.690,.027,.512)
    p$point_size<-.64;p$show_gene_names<-TRUE
    print(p)
  }
  figure2 <- function() {
    page()
    text(.02,.986,"TASSEL maize tutorial: quality control and regional LD",cex=1.34,font=2,adj=c(0,1))
    text(.02,.955,"Official rTASSEL data  |  281 samples  |  10 chromosomes  |  AGPv1 coordinates",cex=.97,adj=c(0,1))
    graphics::segments(.02,.930,.98,.930,col=rule,lwd=.7)
    label(.02,.913,"a","Auditable marker filtering")
    xs<-c(.025,.275,.525,.775);xe<-xs+.202
    titles<-c("3,093","3,045","3,045","2,526")
    details<-c("Input markers","Allele-consistent","Missing fraction ≤0.20","MAF ≥0.05")
    for(i in 1:4) {
      graphics::rect(xs[i],.793,xe[i],.869,col=if(i==4)"#ECF5F3" else "#F7F9FA",border=rule,lwd=.6)
      text((xs[i]+xe[i])/2,.844,titles[i],cex=1.5,font=2)
      text((xs[i]+xe[i])/2,.813,details[i],cex=.97)
      if(i<4)arrow(xe[i]+.009,.832,xs[i+1]-.011,.832)
    }
    text(.125,.775,"48 conflicts excluded",cex=.89,adj=c(.5,1))
    text(.625,.775,"0 removed at missingness step",cex=.89,adj=c(.5,1))
    text(.875,.775,"519 low-MAF markers excluded",cex=.89,adj=c(.5,1))
    label(.02,.734,"b","Retention across chromosomes")
    label(.686,.734,"c","Retained allele frequencies")
    old<-graphics::par(no.readonly=TRUE)
    graphics::par(fig=c(.028,.641,.468,.703),mar=c(2.6,3.8,.3,.3),new=TRUE,cex=.9)
    z<-result$bychr
    graphics::barplot(rbind(z$retained_markers,z$allele_consistent-z$retained_markers,z$raw_markers-z$allele_consistent),
      names.arg=z$chromosome,col=c(teal,"#C8D0D4",orange),border=NA,ylim=c(0,560),las=1,cex.axis=.87,cex.names=.88)
    graphics::mtext("Markers",side=2,line=2.65,cex=.9)
    graphics::mtext("Chromosome",side=1,line=1.65,cex=.9)
    graphics::par(fig=c(.684,.982,.475,.703),mar=c(2.8,3.2,.3,.3),new=TRUE,cex=.9)
    graphics::hist(result$maize$variants$maf,breaks=seq(.05,.5,.05),col=blue,border="white",
      xlab="",ylab="",main="",las=1,cex.axis=.86,xlim=c(.05,.5),axes=FALSE)
    graphics::axis(1,at=c(.05,.2,.35,.5),labels=c("0.05","0.20","0.35","0.50"),cex.axis=.86)
    graphics::axis(2,las=1,cex.axis=.86)
    graphics::mtext("MAF",side=1,line=1.85,cex=.9)
    graphics::mtext("Markers",side=2,line=2.15,cex=.9)
    graphics::par(old)
    graphics::par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=TRUE);graphics::plot.new();graphics::plot.window(c(0,1),c(0,1))
    for(i in 1:3) {
      x<-c(.079,.273,.434)[i]
      graphics::rect(x,.451,x+.012,.461,col=c(teal,"#C8D0D4",orange)[i],border=NA)
      text(x+.019,.456,c("Retained","Low MAF","Allele conflict")[i],cex=.9,adj=0)
    }
    text(.708,.456,"Median MAF = 0.2528",cex=.95,adj=0)
    graphics::segments(.02,.433,.98,.433,col=rule,lwd=.7)
    label(.02,.418,"d","Objectively selected 1-Mb region")
    text(.058,.386,"Chr 1:267,615,478–268,615,477  |  16 markers; 120 pairs; pairwise n = 226–277",cex=.98,adj=c(0,1))
    p<-result$maize_plot;p$title<-NULL;p$viewport<-c(.02,.98,.030,.371)
    p$position_scale<-"index";p$grid_width<-.4
    print(p)
    graphics::par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=TRUE);graphics::plot.new();graphics::plot.window(c(0,1),c(0,1))
    text(.979,.023,"Equal SNP spacing; dosage r²; median = 0.0353",cex=.91,adj=c(1,0))
  }
  draw_file <- function(stem,draw,width,height) {
    for(ext in c("pdf","svg","png")) {
      path<-file.path(output,paste0(stem,".",ext))
      if(ext=="pdf")grDevices::cairo_pdf(path,width=width,height=height,pointsize=8,family="sans")
      else if(ext=="svg")grDevices::svg(path,width=width,height=height,pointsize=8,family="sans")
      else grDevices::png(path,width=width,height=height,units="in",res=600,pointsize=8,type="cairo")
      tryCatch(draw(),finally=grDevices::dev.off())
    }
  }
  draw_file("Figure_1_workflow",figure1,7.5,6.7)
  draw_file("Figure_2_reference_results",figure2,7.5,7.0)
  invisible(output)
}
