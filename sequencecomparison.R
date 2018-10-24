#!/usr/bin/R

library(seqinr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Requires three arguments: 1. a Rscript generated fasta file, 2. an mpileup generated fasta file, and 3. a BTV8 reference fasta file. Please retry providing these arguments.")
}

sequencesR<-read.fasta(args[1])
sequencesM<-read.fasta(args[2])
sequencesRef<-read.fasta(args[3])

workingnames<-names(sequencesR)

k<-1
difference<-vector("list")
for (i in 1:length(sequencesR)) {
  for (j in 1:length(sequencesR[[i]])) {
    if (sequencesR[[which(names(sequencesR)==names(sequencesR)[i])]][j]!=sequencesM[[which(names(sequencesM)==names(sequencesR)[i])]][j]) {
      workingdifference<-c(as.character(names(sequencesR)[i]),as.character(j),as.character(sequencesR[[which(names(sequencesR)==names(sequencesR)[i])]][j]),as.character(sequencesM[[which(names(sequencesM)==names(sequencesR)[i])]][j]),as.character(sequencesRef[[which(names(sequencesRef)==names(sequencesR)[i])]][j]))
      difference[[k]]<-workingdifference
      k<-k+1
    }
  }
}

if (!is.null(difference[1][[1]])) {
  message("Different")
  differences<-do.call("rbind",lapply(difference, function(x) {return(x)}))
  colnames(differences)<-c("Segment","Position","RscriptConsensus","mpileupConsensus","Reference")
  write.csv(differences,file="differences.csv",row.names=F)
} else {
  message("Same.")
}
