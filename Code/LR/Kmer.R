#**************************
# Transfer sequence to k-mer counts
# Tian Tian
# tt72@njit.edu
#**************************

rm(list=ls())
library(Biostrings)
library(rhdf5)
library(foreach)
library(doMC)
registerDoMC()

#**************************
# load sequence
#**************************
len = 10
index = 1

train.pos.fa = readDNAStringSet(paste("m6A_", len, "_train.pos_", index, ".fa", sep=""))
test.pos.fa = readDNAStringSet(paste("m6A_", len, "_test.pos_", index, ".fa", sep=""))
valid.pos.fa = readDNAStringSet(paste("m6A_", len, "_valid.pos_", index, ".fa", sep=""))

train.neg.fa = readDNAStringSet(paste("m6A_", len, "_train.neg_", index, ".fa", sep=""))
test.neg.fa = readDNAStringSet(paste("m6A_", len, "_test.neg_", index, ".fa", sep=""))
valid.neg.fa = readDNAStringSet(paste("m6A_", len, "_valid.neg_", index, ".fa", sep=""))

#**************************
# calculate the frequency of oligos with length of 1-6bp
#**************************
train.pos.1mer = oligonucleotideFrequency(train.pos.fa, 1)
test.pos.1mer = oligonucleotideFrequency(test.pos.fa, 1)
valid.pos.1mer = oligonucleotideFrequency(valid.pos.fa, 1)

train.pos.2mer = oligonucleotideFrequency(train.pos.fa, 2)
test.pos.2mer = oligonucleotideFrequency(test.pos.fa, 2)
valid.pos.2mer = oligonucleotideFrequency(valid.pos.fa, 2)

train.pos.3mer = oligonucleotideFrequency(train.pos.fa, 3)
test.pos.3mer = oligonucleotideFrequency(test.pos.fa, 3)
valid.pos.3mer = oligonucleotideFrequency(valid.pos.fa, 3)

train.pos.4mer = oligonucleotideFrequency(train.pos.fa, 4)
test.pos.4mer = oligonucleotideFrequency(test.pos.fa, 4)
valid.pos.4mer = oligonucleotideFrequency(valid.pos.fa, 4)

train.pos.5mer = oligonucleotideFrequency(train.pos.fa, 5)
test.pos.5mer = oligonucleotideFrequency(test.pos.fa, 5)
valid.pos.5mer = oligonucleotideFrequency(valid.pos.fa, 5)

train.pos.6mer = oligonucleotideFrequency(train.pos.fa, 6)
test.pos.6mer = oligonucleotideFrequency(test.pos.fa, 6)
valid.pos.6mer = oligonucleotideFrequency(valid.pos.fa, 6)

train.pos.kmer = cbind(train.pos.1mer, train.pos.2mer, train.pos.3mer, train.pos.4mer, train.pos.5mer, train.pos.6mer)
test.pos.kmer = cbind(test.pos.1mer, test.pos.2mer, test.pos.3mer, test.pos.4mer, test.pos.5mer, test.pos.6mer)
valid.pos.kmer = cbind(valid.pos.1mer, valid.pos.2mer, valid.pos.3mer, valid.pos.4mer, valid.pos.5mer, valid.pos.6mer)

train.neg.1mer = oligonucleotideFrequency(train.neg.fa, 1)
test.neg.1mer = oligonucleotideFrequency(test.neg.fa, 1)
valid.neg.1mer = oligonucleotideFrequency(valid.neg.fa, 1)

train.neg.2mer = oligonucleotideFrequency(train.neg.fa, 2)
test.neg.2mer = oligonucleotideFrequency(test.neg.fa, 2)
valid.neg.2mer = oligonucleotideFrequency(valid.neg.fa, 2)

train.neg.3mer = oligonucleotideFrequency(train.neg.fa, 3)
test.neg.3mer = oligonucleotideFrequency(test.neg.fa, 3)
valid.neg.3mer = oligonucleotideFrequency(valid.neg.fa, 3)

train.neg.4mer = oligonucleotideFrequency(train.neg.fa, 4)
test.neg.4mer = oligonucleotideFrequency(test.neg.fa, 4)
valid.neg.4mer = oligonucleotideFrequency(valid.neg.fa, 4)

train.neg.5mer = oligonucleotideFrequency(train.neg.fa, 5)
test.neg.5mer = oligonucleotideFrequency(test.neg.fa, 5)
valid.neg.5mer = oligonucleotideFrequency(valid.neg.fa, 5)

train.neg.6mer = oligonucleotideFrequency(train.neg.fa, 6)
test.neg.6mer = oligonucleotideFrequency(test.neg.fa, 6)
valid.neg.6mer = oligonucleotideFrequency(valid.neg.fa, 6)

train.neg.kmer = cbind(train.neg.1mer, train.neg.2mer, train.neg.3mer, train.neg.4mer, train.neg.5mer, train.neg.6mer)
test.neg.kmer = cbind(test.neg.1mer, test.neg.2mer, test.neg.3mer, test.neg.4mer, test.neg.5mer, test.neg.6mer)
valid.neg.kmer = cbind(valid.neg.1mer, valid.neg.2mer, valid.neg.3mer, valid.neg.4mer, valid.neg.5mer, valid.neg.6mer)

train.kmer = rbind(train.pos.kmer, train.neg.kmer)
test.kmer = rbind(test.pos.kmer, test.neg.kmer)
valid.kmer = rbind(valid.pos.kmer, valid.neg.kmer)

train.y = c(rep(1, nrow(train.pos.kmer)), rep(0, nrow(train.neg.kmer)))
test.y = c(rep(1, nrow(test.pos.kmer)), rep(0, nrow(test.neg.kmer)))
valid.y = c(rep(1, nrow(valid.pos.kmer)), rep(0, nrow(valid.neg.kmer)))

#**************************
# write k-mer counts to file in hdf5 format
#**************************
h5createFile(paste("LR_m6A_", len, ".train_", index, ".hdf5", sep=""))
h5createFile(paste("LR_m6A_", len, ".valid_", index, ".hdf5", sep=""))
h5createFile(paste("LR_m6A_", len, ".test_", index, ".hdf5", sep=""))

h5createDataset(paste("LR_m6A_", len, ".train_", index, ".hdf5", sep=""), "X", c(nrow(train.kmer), ncol(train.kmer)),
                storage.mode="integer", chunk=c(1, ncol(train.kmer)), level=9)
h5write(train.kmer, paste("LR_m6A_", len, ".train_", index, ".hdf5", sep=""), "X")
h5write(train.y, paste("LR_m6A_", len, ".train_", index, ".hdf5", sep=""), "Y")
h5write(valid.kmer, paste("LR_m6A_", len, ".valid_", index, ".hdf5", sep=""), "X")
h5write(valid.y, paste("LR_m6A_", len, ".valid_", index, ".hdf5", sep=""), "Y")
h5write(test.kmer, paste("LR_m6A_", len, ".test_", index, ".hdf5", sep=""), "X")
h5write(test.y, paste("LR_m6A_", len, ".test_", index, ".hdf5", sep=""), "Y")
