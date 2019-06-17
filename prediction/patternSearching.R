#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=2){ 
 stop("incorrect number of args")
}

library(Biostrings)

x <- function(file1 = args[1], out = args[2]){
  fasta <- readDNAStringSet(file1)
  fileConn <- file(out, "w")
  
  seq <- "T[AC]A..[AG]TGA[CT]...GC[AG][AT]{4}"
  z <- grep(seq, fasta)

  write(c(length(fasta), names(fasta[z])), fileConn)
  close(fileConn)
}

x()
