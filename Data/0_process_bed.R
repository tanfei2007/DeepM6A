rm(list=ls())
library(Biostrings)
setwd("~/Documents/Research/m6A/Drosophila/dm6/")

seq.dr = readDNAStringSet("./dm6_genome/dm6.fa")

seqLength = lapply(seq.dr, function(z){nchar(toString(z))})
seqLength.d = as.data.frame(seqLength)
seqLength.d = t(seqLength.d)

m6A = read.table("825_m6A_25x_novar_dev03.bed", sep="\t", stringsAsFactors=F)
m6A.bed = data.frame(chrom=m6A$V1, chromStart=m6A$V2, chromEnd=m6A$V3, 
                     name=".", score=".", strand=m6A$V4)
m6A.bed$chromStart = as.integer(apply(m6A, 1, function(z){
  if(z[1]=="chrM") {
    startPos = as.numeric(z[2])-1010
    if(startPos>=0)
      return(startPos)
    else
      return(0)
  } else {
    startPos = as.numeric(z[2])-1010
    if(startPos>=0)
      return(startPos)
    else
      return(0)
  }
}))
m6A.bed$chromEnd = as.integer(apply(m6A, 1, function(z){
  if(z[1]=="chrM") {
    endPos = as.numeric(z[2])+1010
    length = seqLength[z[1]]
    if(endPos<=length)
      return(endPos)
    else
      return(length)
  } else {
    endPos = as.numeric(z[2])+1010
    length = seqLength[z[1]]
    if(endPos<=length)
      return(endPos)
    else
      return(length)
  }
}))
write.table(m6A.bed, "6mA_complement.bed", quote=F, col.names=F, row.names=F, sep="\t")
