#!/usr/bin/R

library(stringr)
library(ggplot2)

##generate per segment TPM barplot
tpms<-read.table("results.xprs",header=T)
temp<-do.call("rbind",str_split(as.character(tpms$target_id),"seg"))
tpms$serotype<-as.factor(temp[,1])
tpms$segment<-as.factor(temp[,2])
ggplot(tpms,aes(serotype,tpm)) + geom_bar(stat="identity") + facet_grid(rows = vars(segment), scales = "free") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("BTV Serotype") + ylab("TPM")
ggsave("figuresBTV.pdf")
