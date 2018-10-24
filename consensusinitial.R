#!/usr/bin/R
library(vcfR)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Requires two arguments: 1. a vcf file and 2. an associated fasta reference file. Please retry providing these arguments.")
}

##read in vcf and reference sequences
vcf<-read.vcfR(args[1])
sequences<-read.fasta(args[2])

if (nrow(vcf@fix)>0) {
  ##extract data from vcf
  data<-as.data.frame(vcf@fix)
  
  ##filter for variants that passed verification
  data<-data[as.character(data$FILTER)%in%"PASS",]
  data$INFO<-as.character(data$INFO)
  data$POS<-as.numeric(as.character(data$POS))
  
  ##extract relevant metadata (allele freq plus depth at variant)
  tempINFOmat<-do.call("rbind",strsplit(data$INFO,";"))
  data$Depth<-as.numeric(do.call("rbind",strsplit(tempINFOmat[,1],"="))[,2])
  data$AF<-as.numeric(do.call("rbind",strsplit(tempINFOmat[,2],"="))[,2])
  data$chrpos<-paste(data$CHROM,data$POS,sep="_")
  
  ##generate list of unique sites
  variablepositions<-vector("list",length(unique(data$chrpos)))
  
  ##for pre-final base calls replace with dominant allele
  for (i in 1:length(unique(data$chrpos))) {
    subdata<-data[data$chrpos%in%unique(data$chrpos)[i],]
    workinglist<-vector("list",length=3)
    workinglist[[1]]<-subdata$CHROM[1]
    workinglist[[2]]<-subdata$POS[1]
    if (sum(subdata$AF)<0.5) {
      workinglist[[3]]<-as.character(subdata$REF[1])
    } else {
      ##if non reference allele dominant call non reference allele
      if (((1-sum(subdata$AF))>subdata$AF[which.max(subdata$AF)])==1) {
        workinglist[[3]]<-as.character(subdata$REF[1])
      } else {
	workinglist[[3]]<-as.character(subdata$ALT[which.max(subdata$AF)])
      }
    }
    variablepositions[[i]]<-workinglist
  }
  
  ##replace in sequences based on correct mapping
  
  for (i in 1:length(variablepositions)) {
    base<-NA
    base<-variablepositions[[i]][[3]]
    sequences[[which(names(sequences)==as.character(variablepositions[[i]][[1]]))]][variablepositions[[i]][[2]]]<-base
  }
  
  ##write sequences out
  
  write.fasta(sequences,names(sequences),args[2])
}
