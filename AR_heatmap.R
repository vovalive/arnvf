rm(list=ls(all=T))

library(ggplot2)
library(gplots)
library(RColorBrewer)

coeff=1000

AR.counts.all <- read.delim("AR.counts.all.txt", stringsAsFactors=FALSE)
AR.covered.100pr <- read.delim("AR.covered.100pr.txt", header=FALSE, stringsAsFactors=FALSE)
AR.genes.desc <- read.delim2("AR.genes.desc.txt", stringsAsFactors=FALSE,quote="")


AR=AR.counts.all
row.names(AR)=AR$name
AR$name=NULL
AR=data.matrix(AR)
colnames(AR)
head(AR)

AR=AR[,-1]

for(i in 1:ncol(AR)){
  AR[,i]=AR[,i]/sum(AR[,i])*coeff
}
colnames(AR)
colnames(AR)=gsub(x = colnames(AR),pattern = ".fastq.AR.sorted.bam",replacement = "")

AR.100=AR[which(row.names(AR) %in% AR.covered.100pr$V2),]
plot(hclust(dist(t(AR.100))))
heatmap.2(log2(AR.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,10))

a=heatmap.2(log2(AR.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,10))
str(a)
genelist=as.data.frame(rev(colnames(a$carpet)))
names(genelist)='name'
description=merge(genelist,AR.genes.desc,by='name',all.x = T,sort=F)

pdf(file = 'AR.heatmap.pdf',width = 30,height = 15)
heatmap.2(log2(AR.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,140),labRow = description$description,keysize = .5,cexRow = 1.1)
dev.off()
