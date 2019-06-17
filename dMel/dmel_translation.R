#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=2){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)
library(seqinr)

#https://bioinfo.umassmed.edu/bootstrappers/guides/main/r_writeFasta.html
writeFasta<-function(data, fileConn){
  fastaLines = c() 
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }

  write(fastaLines, fileConn, append = TRUE)
}


x <- function(file1 = args[1], out = args[2]){
 fasta <- read.fasta(file1)
 fileConn <- file(out, "w")

 n <- names(fasta)
 v1 <- c()
 v2 <- c()

 for(i in 1:length(n)){
    v1[i] <- n[i]
    v2[i] <- toupper(paste(fasta[[i]], collapse = ""))
 }
  
 data = dplyr::data_frame(name = v1, seq = v2)
 writeFasta(data, fileConn)

 close(fileConn)

}

x()

