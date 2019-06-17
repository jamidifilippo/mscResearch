#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=2){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)
library(seqinr)
library(stringr)

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

  n <- word(names(fasta))
  v1 <- c()
  v2 <- c()

  j <- 1

  for(i in 1:length(n)){
    if(subseq(n[i], end = nchar(n[i]), width = 3) == ".t1"){
      v1[j] <- n[i]
    
      #remove trailing '*'
      x <- translate(fasta[[i]])
      v2[j] <- subseq(paste(x, collapse = ""), start=1, width=length(x)-1)
    
      j <- j+1
    }
  }
  
  data = dplyr::data_frame(name = v1, seq = v2)
  writeFasta(data, fileConn)

  close(fileConn)
}

x()

