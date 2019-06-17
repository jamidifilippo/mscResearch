#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=3){ 
 stop("incorrect number of args")
}

library(rtracklayer)
library(Biostrings)
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


x <- function(file1 = args[1], feature=args[2], out=args[3]){
  #open file conn
  fileConn <- file(out, "w")

  #read in gff(3) files
  files <- list.files(path = "hydra2.0_genemodels.gff3/", full.names = TRUE)
  numGFF <- length(files)
  fa <- readDNAStringSet(file1)

  for(i in 1:numGFF){
    gff <- readGFF(files[i])

    #get all 'genes' in this file & remove misc cols
    genes <- gff[which(gff$type==feature),c(1,3,4,5,7,10)]
    #get fasta seq corresponding to the genes in this gff(3)
    fasta <- fa[genes$seqid[1]]

    v1 <- genes$ID
    v2 <- c()
    v3 <- c()
  
    start <- genes$start
    stop <- genes$end
  
    for(j in 1:nrow(genes)){
    
      if(genes$strand[j]=='+'){
        if(start[j]-1000<=0){
          v2[j] <- toString(subseq(fasta, end = start[j]-1))
	  v3[j] <- toString(c("1", start[j]-1))
        }else{
          v2[j] <- toString(subseq(fasta, start = start[j]-1000, width = 1000))
	  v3[j] <- toString(c(start[j]-1000, start[j]-1))
        }
      
      }else{
        if(stop[j]+1000 >= width(fasta)){
          x <- subseq(fasta, start = stop[j]+1)
	  v3[j] <- toString(c(stop[j]+1, width(fasta)))
        }else{
          x <- subseq(fasta, end = stop[j]+1000, width = 1000)
	  v3[j] <- toString(c(stop[j]+1000, stop[j]+1))
        }
   
        v2[j] <- toString(reverseComplement(x))
      }
    }
  
    n <- c()

    for(k in 1: length(v1)){
      n <- c(n, paste(c(v1[k], v3[k]) , collapse = " "))
    }


    require(dplyr)
  
    data = dplyr::data_frame(name = n, seq = v2)
    writeFasta(data, fileConn)
  }

 close(fileConn)
}

x()

