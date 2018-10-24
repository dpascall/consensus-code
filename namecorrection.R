#!/usr/bin/R
library(seqinr)
args = commandArgs(trailingOnly=TRUE)

sequences<-read.fasta(args[1])

snames<-paste(args[2],names(sequences),sep="_")

write.fasta(sequences, snames, args[1])
