rm(list=ls(all=T))

library(ggplot2)
library(gplots)
library(RColorBrewer)

coeff=1000

VF.counts.all <- read.delim("VF.counts.all.txt", stringsAsFactors=FALSE)
VF.covered.100pr <- read.delim("VF.covered.100pr.txt", header=FALSE, stringsAsFactors=FALSE)
VF.genes.desc <- read.delim2("VF.genes.desc.txt", stringsAsFactors=FALSE,quote="")


VF=VF.counts.all
row.names(VF)=VF$name
VF$name=NULL
VF=data.matrix(VF)
colnames(VF)
head(VF)

VF=VF[,-1]

for(i in 1:ncol(VF)){
  VF[,i]=VF[,i]/sum(VF[,i])*coeff
}
colnames(VF)
colnames(VF)=gsub(x = colnames(VF),pattern = ".fastq.VF.sorted.bam",replacement = "")

VF.100=VF[which(row.names(VF) %in% VF.covered.100pr$V2),]
plot(hclust(dist(t(VF.100))))
heatmap.2(log2(VF.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,10))

a=heatmap.2(log2(VF.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,10))
str(a)
genelist=as.data.frame(rev(colnames(a$carpet)))
names(genelist)='name'
description=merge(genelist,VF.genes.desc,by='name',all.x = T,sort=F)

pdf(file = 'VF.heatmap.pdf',width = 30,height = 15)
heatmap.2(log2(VF.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(8,90),labRow = description$description,keysize = .5,cexRow = 1.3,cexCol = 3)
dev.off()
