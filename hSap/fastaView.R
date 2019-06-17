#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=1){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)

x <- function(file=args[1]){
 fasta <- readDNAStringSet(file)
 print(head(fasta))
 print(length(fasta))
 print(names(fasta)[1])
}

x()
