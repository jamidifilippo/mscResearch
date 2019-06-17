#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=1){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)

x <- function(file=args[1]){
 gff <- readGFF(file)
 print(head(gff))
}

x()
