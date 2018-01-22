rm(list=ls())
setwd("~/Documents/Research/m6A/Drosophila/dm6/")
library(Biostrings)
library(stringi)
library(foreach)
library(doMC)
registerDoMC()

seq.dr = readDNAStringSet("./dm6_genome/dm6.fa")

seqLength = lapply(seq.dr, function(z){nchar(toString(z))})
seqLength.d = as.data.frame(seqLength)
seqLength.d = t(seqLength.d)
write.table(seqLength.d, "dm6.chrom.sizes.txt", col.names=F, sep="\t", quote=F)

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

m6A = read.table("825_m6A_25x_novar_dev03.bed", sep="\t", stringsAsFactors=F)

m6A.complement.bed = read.table("6mA.complement.bed", sep="\t", quote="", stringsAsFactors=F)
negative.seq = data.frame(chr=rep(NA, nrow(m6A)), coordinate=NA, strand=NA)
seqLength.d1 = as.data.frame(seqLength.d[row.names(seqLength.d) %in% unique(m6A$V1),,drop=F])
seqLength.d1$V1 = as.numeric(seqLength.d1$V1)
seqLength.d1$V2 = seqLength.d1$V1/sum(seqLength.d1$V1)

chrs = rep(rownames(seqLength.d1), times=nrow(m6A)*seqLength.d1$V2)
chrs = c(chrs, sample(unique(rownames(seqLength.d1)), nrow(negative.seq)-length(chrs), replace=T))
chrs = sort(chrs)
negative.seq$chr = chrs

chr.table = table(negative.seq$chr)
coors = c()
strand = c()
for(i in 1:length(chr.table)) {
  chr = names(chr.table)[i]
  m6A.complement.bed.chr = m6A.complement.bed[m6A.complement.bed$V1==chr,]
  seqs = unlist(mapply(seq, m6A.complement.bed.chr[,2], m6A.complement.bed.chr[,3]))
  seq.chr = c()
  number = as.numeric(chr.table[i])
  number0 = 0
  strand.chr = c()
  while(T) {
    seq.chr.try = sample(seqs, number)
    characters = unlist(lapply(seq.chr.try, function(z) {toString(subseq(seq.dr[chr], z, z))}))
    characters.t = table(characters)
    if(!is.na(characters.t['A']) & !is.na(characters.t['T'])) {
      number.c = as.numeric(characters.t['A'] + characters.t['T'])
    } else if(!is.na(characters.t['A']) & is.na(characters.t['T'])) {
      number.c = as.numeric(characters.t['A'])
    } else if(is.na(characters.t['A']) & !is.na(characters.t['T'])) {
      number.c = as.numeric(characters.t['T'])
    } else {
      number.c = NA
    }
    if(is.na(number.c))
      next
    if(number0+number.c > number) {
      seq.chr1 = c(seq.chr.try[which(characters=='A')], seq.chr.try[which(characters=='T')])
      seq.chr2 = sample(seq.chr1, number-number0)
      characters2 = unlist(lapply(seq.chr2, function(z) {toString(subseq(seq.dr[chr], z, z))}))
      seq.chr = c(seq.chr, seq.chr2)
      strand.chr1 = rep('+', number-number0)
      strand.chr1[which(characters2=='T')] = '-'
      strand.chr = c(strand.chr, strand.chr1)
      break
    }
    else if(number0+number.c == number) {
      seq.chr = c(seq.chr, seq.chr.try[which(characters=='A')], seq.chr.try[which(characters=='T')])
      if(is.na(as.numeric(characters.t['T']))) {
        strand.chr = c(strand.chr, rep('+', as.numeric(characters.t['A'])))
      } else if(is.na(as.numeric(characters.t['A']))) {
        strand.chr = c(strand.chr, rep('-', as.numeric(characters.t['T'])))
      } else {
        strand.chr = c(strand.chr, rep('+', as.numeric(characters.t['A'])), rep('-', as.numeric(characters.t['T'])))
      }
      break
    }
    else {
      seq.chr1 = c(seq.chr.try[which(characters=='A')], seq.chr.try[which(characters=='T')])
      seq.chr = c(seq.chr, seq.chr1)
      strand.chr = c(strand.chr, rep('+', as.numeric(characters.t['A'])), rep('-', as.numeric(characters.t['T'])))
      number0 = number0 + number.c
      seqs = setdiff(seqs, seq.chr1)
    }
  }
  coors = c(coors, seq.chr)
  strand = c(strand, strand.chr)
}
negative.seq$coordinate = coors
negative.seq$strand = strand


len = 20

negative.seq$context = apply(negative.seq, 1, function(z){
  startPos = as.numeric(z[2])-len
  endPos = as.numeric(z[2])+len
  length = as.numeric(seqLength[z[1]])
  if(startPos>0 & endPos<=length) {
    if(z[3]=="+")
      return(toString(subseq(seq.dr[z[1]], startPos, endPos)))
    else
      return(rev.comp(toString(subseq(seq.dr[z[1]], startPos, endPos))))
  }
  else {
    if(startPos<=0) {
      if(z[3]=="+")
        return(paste(stri_dup('N', 1-startPos), toString(subseq(seq.dr[z[1]], 1, endPos)), sep=""))
      else
        return(paste(rev.comp(toString(subseq(seq.dr[z[1]], 1, endPos))), stri_dup('N', 1-startPos), sep=""))
    }
    if(endPos>length) {
      if(z[3]=="+")
        return(paste(toString(subseq(seq.dr[z[1]], startPos, length)), stri_dup('N', endPos-length), sep=""))
      else
        return(paste(stri_dup('N', endPos-length), rev.comp(toString(subseq(seq.dr[z[1]], startPos, length))), sep=""))
    }
  }
})

setwd("positive_and_negative_sequence/")
m6A.m = data.frame(chr=m6A$V1, coordinate=(m6A$V2+m6A$V3+1)/2, strand=m6A$V4, context=m6A$V5)
write.csv(m6A.m, paste("m6A_",len,".pos.csv",sep=""), row.names=F)
write.csv(negative.seq, paste("m6A_",len,".neg.csv",sep=""), row.names=F)
