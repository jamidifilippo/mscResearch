#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=3){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)

#https://bioinfo.umassmed.edu/bootstrappers/guides/main/r_writeFasta.html
writeFasta<-function(data, fileConn){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }

  write(fastaLines, fileConn, append = TRUE)
}

x <- function(file1 = args[1], file2 = args[2], out = args[3]){
 fasta <- readDNAStringSet(file1)
 fasta2 <- readDNAStringSet(file2)

 fileConn <- file(out, "w")

 v1 <- c()
 v2 <- c()
 v3 <- c()
  
for(i in 1:length(fasta)){

  neg <- gregexpr("complement", names(fasta)[i])[[1]][1]

  #protein name 
  index <- gregexpr(" ", names(fasta)[i])[[1]][1]
  v1[i] <- subseq(names(fasta)[i], start=1, end=index-1)

  #determine +/- strand
  index1 <- gregexpr(":", names(fasta)[i])[[1]][1]  
  value <- subseq(names(fasta)[i], start=index1+1, width=4)

  #determine chrom
  index3 <- gregexpr("loc=", names(fasta)[i])[[1]][1]
  value2 <- subseq(names(fasta)[i], start=index3+4, end=index1)
  
  fa <- fasta2[grep(value2, names(fasta2))][1]
  
  
  if(neg!='-1'){ # - strand
      index3 <- gregexpr("\\.", names(fasta)[i])[[1]][1]
      index4 <- gregexpr("[\\,|\\)]", names(fasta)[i])[[1]][1]
      stop <- subseq(names(fasta)[i], start=index3+2, end=index4-1)
      stop <- as.numeric(stop)
      
      if(stop+1000 >= width(fa)){
        x <- subseq(fa, start = stop+1)
	v3[i] <- toString(c(stop+1, width(fa)))
      }else{
        x <- subseq(fa, end = stop+1000, width = 1000)
	v3[i] <- toString(c(stop+1000, stop+1))
      }   
          
      v2[i] <- toString(reverseComplement(x))
   
  }else{
    index3 <- gregexpr("\\.", names(fasta)[i])[[1]][1]

    if(value=="join"){  # + strand
      index2 <- gregexpr("\\(", names(fasta)[i])[[1]][1]
      start <- subseq(names(fasta)[i], start=index2+1, end=index3-1)
      start <- as.numeric(start)
    }else{
      index4 <- gregexpr(":", names(fasta)[i])[[1]][1]
      start <- subseq(names(fasta)[i], start=index4+1, end=index3-1)
      start <- as.numeric(start)
    }

      if(start-1000<=0){
        v2[i] <- toString(subseq(fa, end = start-1))
	v3[i] <- toString(c("1", start-1))
      }else{
        v2[i] <- toString(subseq(fa, start = start-1000, width = 1000))
	v3[i] <- toString(c(start-1000, start-1))
      }     
    }
}

 n <- c()

 for(k in 1: length(v1)){
   n <- c(n, paste(c(v1[k], v3[k]) , collapse = " "))
 }

 data = dplyr::data_frame(name = n, seq = v2)
 writeFasta(data, fileConn)

 close(fileConn)
}

x()

