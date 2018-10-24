#!/usr/bin/R
library(seqinr)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Requires two arguments: 1. a depth tab separated txt file and 2. an associated fasta reference file. Please retry providing these arguments.")
}

##read in depth and reference sequences
depth<-read.table(args[1],header=F,sep="\t")
depth[,1]<-as.character(depth[,1])
depth[,2]<-as.numeric(as.character(depth[,2]))
depth[,3]<-as.numeric(as.character(depth[,3]))
sequences<-read.fasta(args[2])

##filter for low depth
depth<-depth[depth[,3]<5,]

##replace in sequences based on correct mapping
if (nrow(depth)>0) {
  for (i in 1:nrow(depth)) {
    sequences[[which(names(sequences)==depth[i,1])]][depth[i,2]]<-"N"
  }
}

##write sequences out
write.fasta(sequences,names(sequences),args[2])