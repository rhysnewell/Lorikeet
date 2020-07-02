#!/usr/bin/env Rscript
source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")


args <- commandArgs(trailingOnly = T)
snps <- read.table(paste0(args[1], "_snp_locations.tsv"))
window_size <- as.numeric(args[2])*1000
CMplot(snps,type="p",plot.type="d",bin.size=window_size,chr.den.col=c("darkgreen", "yellow", "red"),file="pdf",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)