#!/usr/bin/env Rscript
args = commandArgs(TRUE)

if(length(args)!=3){ 
 stop("incorrect number of args")
}

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


x <- function(file1=args[1], file2=args[2], out=args[3]){
  fasta <- read.fasta(file1)
  names <- read.delim(file2, header=FALSE)
  fileConn <- file(out, "w")

  y <- names[2:nrow(names),1]
  z <- c()

  for(i in 1:(nrow(names)-1)){
    z <- c(z, strsplit(toString(y[[i]]), " ")[[1]][1])
  }

  index <- c()
  f <- c()

  for(k in 1:length(z)){
    index <- grep(z[k], names(fasta), value=FALSE)  
    s <- toupper(paste(fasta[[index]], collapse = ""))
    f <- c(f, toString(s))
   }

  #fasta <- fasta[index]

  require(dplyr)
  
  data = dplyr::data_frame(name = z, seq = f)
  writeFasta(data, fileConn)
}

x()
