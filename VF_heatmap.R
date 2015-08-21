rm(list=ls(all=T))

library(ggplot2)
library(gplots)
library(RColorBrewer)

coeff=1000

VF.counts.all <- read.delim("VF.counts.all.txt", stringsAsFactors=FALSE)
VF.covered.100pr <- read.delim("VF.covered.100pr.txt", header=FALSE, stringsAsFactors=FALSE)
VF.genes.desc <- read.delim2("VF.genes.desc.txt", stringsAsFactors=FALSE)


vf=VF.counts.all
row.names(vf)=vf$name
vf$name=NULL
vf=data.matrix(vf)
for(i in 1:ncol(vf)){
  vf[,i]=vf[,i]/sum(vf[,i])*coeff
}
colnames(vf)
colnames(vf)=gsub(x = colnames(vf),pattern = ".fastq.VF.sorted.bam",replacement = "")

vf.100=vf[which(row.names(vf) %in% VF.covered.100pr$V2),]
plot(hclust(dist(t(vf.100))))
heatmap.2(log2(vf.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')))

a=heatmap.2(log2(vf.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')))
str(a)
genelist=as.data.frame(rev(colnames(a$carpet)))
names(genelist)='name'
description=merge(genelist,VF.genes.desc,by='name',all.x = T,sort = F)



