rm(list=ls())
setwd("~/Documents/Research/m6A/Drosophila/dm6/")
library(Biostrings)
library(stringi)
library(foreach)
library(doMC)
registerDoMC()

rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}

genome.seq = readDNAStringSet("dm6_genome/dm6.fa")
# names(genome.seq) = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC")

seqLength = lapply(genome.seq, function(z){nchar(toString(z))})
seqLength.d = as.data.frame(seqLength)
seqLength.d = t(seqLength.d)

len = 1
m6A = read.csv("dm6_m6A.pos.csv", stringsAsFactors = F, row.names = 1)
m6A$context = apply(m6A, 1, function(z){
  startPos = as.numeric(z[2])-len
  endPos = as.numeric(z[2])+len
  length = as.numeric(seqLength[z[1]])
  if(startPos>0 & endPos<=length) {
    if(z[3]=="+")
      return(toString(subseq(genome.seq[z[1]], startPos, endPos)))
    else
      return(rev.comp(toString(subseq(genome.seq[z[1]], startPos, endPos))))
  }
  else {
    if(startPos<=0) {
      if(z[3]=="+")
        return(paste(stri_dup('N', 1-startPos), toString(subseq(genome.seq[z[1]], 1, endPos)), sep=""))
      else
        return(paste(rev.comp(toString(subseq(genome.seq[z[1]], 1, endPos))), stri_dup('N', 1-startPos), sep=""))
    }
    if(endPos>length) {
      if(z[3]=="+")
        return(paste(toString(subseq(genome.seq[z[1]], startPos, length)), stri_dup('N', endPos-length), sep=""))
      else
        return(paste(stri_dup('N', endPos-length), rev.comp(toString(subseq(genome.seq[z[1]], startPos, length))), sep=""))
    }
  }
})

negative.seq = data.frame(chr=m6A$chr, coordinate=NA, strand=m6A$strand, stringsAsFactors=F)
for(i in 1:nrow(m6A)) {
  currentNeg = negative.seq[negative.seq$chr==negative.seq[i,"chr"] & negative.seq$strand==negative.seq[i,"strand"],]
  currentM6A = m6A[m6A$chr==negative.seq[i,"chr"] & m6A$strand==negative.seq[i,"strand"],]
  context = 1
  while(T) {
    selectCoor = NA
    seqs = c((m6A[i,"coordinate"]-context-1):(m6A[i,"coordinate"]-context), (m6A[i,"coordinate"]+context):(m6A[i,"coordinate"]+context+1))
    seqs = seqs[(!seqs %in% currentNeg$coordinate) & (!seqs %in% currentM6A$coordinate)]
    characters = unlist(lapply(seqs, function(z) {toString(subseq(genome.seq[negative.seq[i,"chr"]], z, z))}))
    if(negative.seq[i,"strand"]=="+" & sum(characters=="A")>1) {
      selectCoor = sample(seqs[which(characters=="A")], 1)
      break
    } else if(negative.seq[i,"strand"]=="+" & sum(characters=="A")==1) {
      selectCoor =seqs[which(characters=="A")]
      break
    } else if(negative.seq[i,"strand"]=="-" & sum(characters=="T")>1) {
      selectCoor = sample(seqs[which(characters=="T")], 1)
      break
    } else if(negative.seq[i,"strand"]=="-" & sum(characters=="T")==1) {
      selectCoor = seqs[which(characters=="T")]
      break
    } else {
      context = context+1
      rm("seqs", "characters")
    }
  }
  # return(selectCoor)
  negative.seq[i, "coordinate"] = selectCoor
  rm("currentNeg", "currentM6A", "seqs", "characters", "selectCoor", "context")
}
# negative.seq$coordinate = unlist(selectCoors)


negative.seq$context = apply(negative.seq, 1, function(z){
  startPos = as.numeric(z[2])-len
  endPos = as.numeric(z[2])+len
  length = as.numeric(seqLength[z[1]])
  if(startPos>0 & endPos<=length) {
    if(z[3]=="+")
      return(toString(subseq(genome.seq[z[1]], startPos, endPos)))
    else
      return(rev.comp(toString(subseq(genome.seq[z[1]], startPos, endPos))))
  }
  else {
    if(startPos<=0) {
      if(z[3]=="+")
        return(paste(stri_dup('N', 1-startPos), toString(subseq(genome.seq[z[1]], 1, endPos)), sep=""))
      else
        return(paste(rev.comp(toString(subseq(genome.seq[z[1]], 1, endPos))), stri_dup('N', 1-startPos), sep=""))
    }
    if(endPos>length) {
      if(z[3]=="+")
        return(paste(toString(subseq(genome.seq[z[1]], startPos, length)), stri_dup('N', endPos-length), sep=""))
      else
        return(paste(stri_dup('N', endPos-length), rev.comp(toString(subseq(genome.seq[z[1]], startPos, length))), sep=""))
    }
  }
})

setwd("./positive_and_negative_sequence_neighbour/")
write.csv(m6A, paste("m6A_",len,".csv",sep=""), row.names=F)
write.csv(negative.seq, paste("negative_",len,".csv",sep=""), row.names=F)
