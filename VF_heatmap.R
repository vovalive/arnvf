rm(list=ls(all=T))

library(ggplot2)
library(gplots)
library(RColorBrewer)

coeff=1000

VF.counts.all <- read.delim("VF.counts.all.txt", stringsAsFactors=FALSE)
VF.covered.100pr <- read.delim("VF.covered.100pr.txt", header=FALSE, stringsAsFactors=FALSE)
VF.genes.desc <- read.delim2("VF.genes.desc.txt", stringsAsFactors=FALSE,quote="")


VF=VF.counts.all
VF=VF[-nrow(VF),]
row.names(VF)=VF.genes.desc$description
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

VF.100=VF[which(VF.counts.all$name %in% VF.covered.100pr$V2),]
plot(hclust(dist(t(VF.100))))
heatmap.2(log2(VF.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(5,10))

pdf(file = 'VF.heatmap.pdf',width = 30,height = 15)
heatmap.2(log2(VF.100+1),trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(8,90),keysize = .5,cexRow = 1.3,cexCol = 3)
dev.off()

### now 0 and 1
VF.100.01=VF.100
VF.100.01[which(VF.100.01 > 0)]=1


pdf(file = 'VF.heatmap.01.pdf',width = 30,height = 15)
heatmap.2(VF.100.01,trace='none',col=rev(brewer.pal(9,'Spectral')),margins=c(8,90),keysize = .5,cexRow = 1.3,cexCol = 3)
dev.off()

write.table(VF.100.01,file = 'VF.covered100percent.txt',sep='\t',quote=F,col.names = NA)