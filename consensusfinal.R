#!/usr/bin/R
library(vcfR)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Requires two arguments: 1. a vcf file and 2. an associated fasta reference file. Please retry providing these arguments.")
}

##create mapping between bases present and IUPAC ambiguity codes - function found on StackOverflow https://stackoverflow.com/questions/27912800/check-whether-two-vectors-contain-the-same-unordered-elements-in-r


IUPAC<-vector("list",length=11)
IUPAC[[1]]<-list("R",c("A","G"))
IUPAC[[2]]<-list("Y",c("C","T"))
IUPAC[[3]]<-list("S",c("C","G"))
IUPAC[[4]]<-list("W",c("A","T"))
IUPAC[[5]]<-list("K",c("G","T"))
IUPAC[[6]]<-list("M",c("A","C"))
IUPAC[[7]]<-list("B",c("C","G","T"))
IUPAC[[8]]<-list("D",c("A","G","T"))
IUPAC[[9]]<-list("H",c("A","C","T"))
IUPAC[[10]]<-list("V",c("A","C","G"))
IUPAC[[11]]<-list("N",c("A","C","G","T"))
IUPACref<-c("R","Y","S","W","K","M","B","D","H","V","N")

SameElements<-function(a, b) return(identical(sort(a), sort(b)))

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
  
  ##for each position, make a list of all bases with frequency greater than 5% at that position (assuming the dominant has a frequency of less than 75%)
  for (i in 1:length(unique(data$chrpos))) {
    subdata<-data[data$chrpos%in%unique(data$chrpos)[i],]
    workinglist<-vector("list",length=3)
    workinglist[[1]]<-subdata$CHROM[1]
    workinglist[[2]]<-subdata$POS[1]
    if (sum(subdata$AF)<0.25) {
      workinglist[[3]]<-as.character(subdata$REF[1])
    } else {
      ##if non reference allele dominant call non reference allele
      if (sum(subdata$AF>=0.75)==1) {
        workinglist[[3]]<-as.character(subdata$ALT[which(subdata$AF>=0.75)])
        ##otherwise make list of all alleles with AF greater than 
      } else {
        subdata<-subdata[subdata$AF>=0.05,]
        tempbases<-rep(NA,nrow(subdata)+1)
	if (nrow(subdata)>0) {
          for (j in 1: nrow(subdata)) {
            tempbases[j]<-as.character(subdata$ALT[j])
          }
	}
        tempbases[length(tempbases)]<-as.character(subdata$REF[1])
        workinglist[[3]]<-tempbases
      }
    }
  variablepositions[[i]]<-workinglist
  }

  ##replace in sequences based on correct mapping
  
  for (i in 1:length(variablepositions)) {
    base<-NA
    if (length(variablepositions[[i]][[3]])==1) {
      base<-variablepositions[[i]][[3]]
    } else {
      for (j in 1:11) {
        if (SameElements(IUPAC[[j]][[2]],variablepositions[[i]][[3]])) {
          base<-IUPAC[[j]][[1]]
          break()
        }
      }
    }
    sequences[[which(names(sequences)==as.character(variablepositions[[i]][[1]]))]][variablepositions[[i]][[2]]]<-base
  }
  
  ##write sequences out
  
  write.fasta(sequences,names(sequences),args[2])
}
