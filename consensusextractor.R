library(seqinr)

##extract directories where consensus sequences stored
folders<-list.dirs()[grep("BTV8final",list.dirs())[!grep("BTV8final",list.dirs())%in%grep("_html_results",list.dirs())]]

##get the names of each sequence stripping "Bld" when present
workingnames<-do.call("rbind",strsplit(folders,"/"))[,2]
names<-do.call("c",lapply((strsplit(workingnames,"_"), function (x) {
	if (x[1]=="Bld") {
		return(x[2])
	} else {
		return(x[1])
	}
})

##get names for reading in of consensus sequences
namesforreadin<-do.call("c",strsplit(workingnames,"_working"))

##read in sequences and add relevant sequences to the right combined file
for (i in 1:length(folders)) {
	if (i==1) {
		seg1<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg2<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg3<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg4<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg5<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg6<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg7<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg8<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg9<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg10<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg1<-seg1[grep(names(seg1),"seg1")[!grep(names(seg1),"seg1")%in%grep(names(seg1),"seg10")]]
		seg2<-seg2[grep(names(seg2),"seg2")]
		seg3<-seg3[grep(names(seg3),"seg3")]
		seg4<-seg4[grep(names(seg4),"seg4")]
		seg5<-seg5[grep(names(seg5),"seg5")]
		seg6<-seg6[grep(names(seg6),"seg6")]
		seg7<-seg7[grep(names(seg7),"seg7")]
		seg8<-seg8[grep(names(seg8),"seg8")]
		seg9<-seg9[grep(names(seg9),"seg9")]
		seg10<-seg10[grep(names(seg10),"seg10")]
	} else {
		workingsequences<-read.fasta(paste(folders[i],"/",namesforreadin[i],"_consensus.fasta",sep=""))
		seg1[[i]]<-workingsequences[grep(names(workingsequences),"seg1")[!grep(names(workingsequences),"seg1")%in%grep(names(workingsequences),"seg10")]]
		seg2[[i]]<-workingsequences[grep(names(workingsequences),"seg2")]
		seg3[[i]]<-workingsequences[grep(names(workingsequences),"seg3")]
		seg4[[i]]<-workingsequences[grep(names(workingsequences),"seg4")]
		seg5[[i]]<-workingsequences[grep(names(workingsequences),"seg5")]
		seg6[[i]]<-workingsequences[grep(names(workingsequences),"seg6")]
		seg7[[i]]<-workingsequences[grep(names(workingsequences),"seg7")]
		seg8[[i]]<-workingsequences[grep(names(workingsequences),"seg8")]
		seg9[[i]]<-workingsequences[grep(names(workingsequences),"seg9")]
		seg10[[i]]<-workingsequences[grep(names(workingsequences),"seg10")]
	}
}

##write files
write.fasta(seg1,names,"~/CombinedSequences/BTV8Segment1.fasta")
write.fasta(seg2,names,"~/CombinedSequences/BTV8Segment2.fasta")
write.fasta(seg3,names,"~/CombinedSequences/BTV8Segment3.fasta")
write.fasta(seg4,names,"~/CombinedSequences/BTV8Segment4.fasta")
write.fasta(seg5,names,"~/CombinedSequences/BTV8Segment5.fasta")
write.fasta(seg6,names,"~/CombinedSequences/BTV8Segment6.fasta")
write.fasta(seg7,names,"~/CombinedSequences/BTV8Segment7.fasta")
write.fasta(seg8,names,"~/CombinedSequences/BTV8Segment8.fasta")
write.fasta(seg9,names,"~/CombinedSequences/BTV8Segment9.fasta")
write.fasta(seg10,names,"~/CombinedSequences/BTV8Segment10.fasta")
