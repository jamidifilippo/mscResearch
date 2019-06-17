#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=4){ 
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


x <- function(file1 = args[1], file2 = args[2], feature=args[3], out=args[4]){
 gff <- readGFF(file1)
 fasta <- readDNAStringSet(file2)

 gff <- gff[which(gff$type==feature),]
 fileConn <- file(out, "w")

 v1 <- gff$Name
 v2 <- c()
 v3 <- c()

 start <- gff$start
 stop <- gff$end
  
 f <- word(names(fasta))
# index <- 1

 for(i in 1:nrow(gff)){
 index <- -1
  
  for(j in 1:length(f)){
	if(gff$seqid[i] == f[j]){
		index <- j
	}
  }

  if(index==-1){next}	

# index <- grep(gff$seqid[i], word(names(fasta)))[1]  
 fa <- fasta[index]

    if(gff$strand[i]=='+'){
      if(start[i]-1000<=0){
        v2[i] <- toString(subseq(fa, end = start[i]-1))
        v3[i] <- toString(c("1", start[i]-1))
      }else{
        v2[i] <- toString(subseq(fa, start = start[i]-1000, width = 1000))
	v3[i] <- toString(c(start[i]-1000, start[i]-1))
      }
      
    }else{
      if(stop[i]+1000 >= width(fa)){
        x <- subseq(fa, start = stop[i]+1)
	v3[i] <- toString(c(stop[i]+1, width(fasta)))
      }else{
        x <- subseq(fa, end = stop[i]+1000, width = 1000)
	v3[i] <- toString(c(stop[i]+1000, stop[i]+1))
      }
    
      v2[i] <- toString(reverseComplement(x))
    }
 }

 n <- c()

 for(k in 1: length(v1)){
   n <- c(n, paste(c(v1[k], v3[k]) , collapse = " "))
 }

 require(dplyr)
 data = dplyr::data_frame(name = n, seq = v2) 
 writeFasta(data, fileConn)

 close(fileConn)
}

x()

