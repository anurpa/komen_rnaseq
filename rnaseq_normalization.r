#Komen project, check FPKM distributions, quantile normalization

#Install packages
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","preprocessCore"))
biocLite("sva")

#Load libraries
library(devtools)
library(Biobase)
library(preprocessCore)
library(dplyr)
library(xlsx)
library(sva)
library(pheatmap)

#Set working directory
setwd("~/Projects/komen/analysis_6jun/")

#Read in batch details
batch=read.xlsx("batch.xlsx",1)
batch=batch[batch$sample %in% colnames(edata3),]

write.table(colnames(edata3),"samples.txt",quote=F,sep="\t")

#Read in expression data
edata=read.table("gene_log_name.txt",sep="\t",header=T)

#Remove samples with low reads mapping<5M, 13 samples
#Atleast 10-15M reads required, for RNA-seq good data
#Most of them were from library BC6M8TACXX with just part of a lane data
colnames(edata)=gsub("_0","",colnames(edata))
edata=edata[,-which(names(edata) %in% c("RNA150213LC_B9","RNA150213LC_H1","RNA150213LC_H2","RNA150213LC_B4","RNA150213LC_D9","RNA150213LC_E4","RNA150213LC_H3","RNA150213LC_D5","RNA150213LC_F5","RNA150213LC_E7","RNA140128LC_184","RNA150213LC_G6","RNA140128LC_54"))]

#remove rows with 0 values
edata=edata[ !rowSums(edata[,colnames(edata)[(2:ncol(edata))]]==0)==ncol(edata)-1, ]

#which genes have duplicate values, aggregate 
dupgene=edata[duplicated(edata[,1]),]
dupdata=edata %>% filter(geneSymbol %in% dupgene$geneSymbol)
dypgene_agg=aggregate(dupdata[,2:104],list(geneSymbol=dupdata$geneSymbol),sum)

edata=edata %>% filter(!geneSymbol %in% dupgene$geneSymbol)
edata=rbind(edata,dypgene_agg)

row.names(edata)=edata[,1]
edata=edata[,-1]

edata1=edata
edata1$zero=apply(edata1,1,function(x) sum(x==0.00000000,na.rm=T))
table(edata1$zero)
#95%non zero values, remove genes with all 0 or no expression data
edata2=edata1[edata1$zero<98,]
edata2=edata2[,-104]
edata3=edata2[,-104]

#Distribution plots
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata2[,1]),col=colramp[1],lwd=3,ylim=c(0,1), main="Expression of all genes",xlab="log2 FPKM")
for(i in 2:20){lines(density(edata2[,i]),lwd=3,col=colramp[i])}

#Quantile normalization
norm_edata = normalize.quantiles(as.matrix(edata3))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.50))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

#Plot after normalization
plot(norm_edata[1,],col=as.factor(batch$batch))

#So the columns of U from the SVD correspond to the principal components x in the PCA. Furthermore, the matrix V from the SVD is equivalent to the rotation matrix returned by prcomp.
svd1 = svd(norm_edata - rowMeans(norm_edata))

#samples vs genes, if genes vs sample use u component, svd1$u
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.factor(batch$batch),pch=16)

combat_edata = ComBat(dat=norm_edata, batch=batch$batch, par.prior=TRUE, prior.plots=FALSE)

svd2= svd(combat_edata - rowMeans(combat_edata))

#samples vs genes, if genes vs sample use u component, svd1$u
plot(svd2$v[,1],svd2$v[,2],xlab="PC1",ylab="PC2",
     col=as.factor(batch$batch),pch=16)

cdata=data.frame(combat_edata)
colnames(cdata)=colnames(edata3)
rownames(cdata)=rownames(edata3)

pamgenes=read.table("pamgenes.txt",header = T)
#Alternate gene names for some of the PAM50 genes in our data, CDCA1=NUF2,KNTC2=ndc80,MIA=CD-RAP,ORC6L=ORC6
pamdata=cdata[rownames(cdata) %in% pamgenes$gene,]
missing=cdata[rownames(cdata) %in% c("NUF2","ORC6","NDC80"),]

cdatawnames=cdata
cdatawnames$gene=rownames(cdata)
mia=filter(cdatawnames,grepl("MIA",rownames(cdatawnames)))

pamdata=rbind(pamdata,missing,mia[,-104])

write.table(pamdata,"pamdata.txt",quote=F,sep="\t")
write.table(cdata,"mRNA_lowcov_removed_quantile_normalized_batch_corrected.txt",quote=F,sep="\t")
