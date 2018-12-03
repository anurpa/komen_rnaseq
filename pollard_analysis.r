##Analysis done for Pollard

#Import libraries
library(dplyr)
library(xlsx)
library(pheatmap)
sessionInfo()

#Set working directory
setwd("~/Projects/komen/analysis_30march2017/")

#Read in expression and annotation data
data1=read.table("mRNA_lowcov_removed_quantile_normalized_batch_corrected.txt",header=T,sep="\t")
annot=read.xlsx("~/Projects/komen/sampletracking/annotation_30march2017.xlsx",1)

#Subset tumors only
tumors=subset(annot, annot$tissue=="Tumor" & annot$confidance!=1)

#kmeans used for csf group classification
csfdata=read.xlsx("csfdata.xlsx",1)
row.names(csfdata)=csfdata[,1]
csfdata=csfdata[,-1]
middle=apply(csfdata,1,median)
csfdata_n=csfdata-middle
csfdata_t=t(csfdata_n)
str(csfdata_t)
set.seed(123)
kc=kmeans(csfdata_t,3,iter=1000)
csf_group=kc$cluster
tumors=cbind(tumors,csf_group)

tumors$csf_exp_group=ifelse(tumors$csf_group==1,"mid",ifelse(tumors$csf_group==2,"low",ifelse(tumors$csf_group==3,"high",NA)))
tumors$csf_exp_group=as.factor(tumors$csf_exp_group)

#Ranges csf k means 
summary(kc$centers[1,])
summary(kc$centers[2,])
summary(kc$centers[3,])

table(tumors$csf_group,tumors$subtype)
annotation_col = data.frame(tumors[,c(2,9)])
rownames(annotation_col)=colnames(csfdata)

#Annotation colors for heatmap
ann_colors = list(
  subtype = c(Basal="red",Normal="darkgreen",LumA="darkblue",LumB="lightblue",Her2="pink"),
 csf_exp_group=c(high="red",mid="orange",low="darkgreen"))

#Check csf group assignment
pheatmap(csfdata,scale = "row",clustering_distance_rows = "correlation",annotation_col = annotation_col)

#TAM markers
markers=read.xlsx("tam_markers_revised.xlsx",1)
tamdata=data1[rownames(data1) %in% markers$tam,colnames(data1) %in% tumors$sample]

annotation_col = data.frame(tumors[,c(2,9)])
rownames(annotation_col)=colnames(tamdata)

subtype_order=order(tumors$subtype)
tamdata_order=tamdata[subtype_order]

setEPS()
postscript("tamdata_heatmap.eps")
pheatmap(tamdata_order,scale = "row",clustering_distance_rows = "correlation",annotation_col = annotation_col,clustering_method = "average",
         annotation_colors = ann_colors,col=colorRampPalette(c('green','black','red'))(299),show_colnames = F,cluster_cols = F)
dev.off()

#median normalize tam data
middle=apply(tamdata,1,median)
tamdata_n=tamdata-middle
tamdata_t=t(tamdata_n)
tamdata_t=as.matrix(tamdata_t)
tamscore=apply(tamdata_t,1,median)

tamscore=as.data.frame(tamscore)
str(tamscore)

tumors$tamscore=tamscore$tamscore

setEPS()
postscript("tam_subtype.eps")
boxplot(tumors$tamscore~tumors$subtype,col=c("red","pink","darkblue","lightblue","darkgreen"),varwidth=T,ylab="TAM signature score",xlab="BRCA subtypes")
text(5,0.4,"p=0.00177")
dev.off()

#ANOVA
summary(aov(tumors$tamscore~tumors$subtype))
aov1=aov(tumors$tamscore~tumors$subtype)
posthoc = TukeyHSD(x=aov1)

#TAM score vs immsig 
summary(aov(tumors$tamscore~tumors$immsig))
aov2=aov(tumors$tamscore~tumors$immsig)
posthoc = TukeyHSD(x=aov2)

#TAM score vs csf
summary(aov(tumors$tamscore~tumors$csf_exp_group))
aov3=aov(tumors$tamscore~tumors$csf_exp_group)
posthoc = TukeyHSD(x=aov3)

setEPS()
postscript("tam_markers_immsig.eps")
par(mfrow=c(2,1))
boxplot(tumors$tamscore~tumors$immsig,col=c("red","pink","darkblue","lightblue"),varwidth=T,ylab="TAM signature score",xlab="Immune signature")
text(4,0.85,"p=1.23e-05")
boxplot(tumors$tamscore~tumors$csf_exp_group,col=c("red","pink","darkblue","lightblue","darkgreen"),varwidth=T,ylab="TAM signature score",xlab="CSF1 signature")
text(3,0.85,"p=7.58e-11")
dev.off()


#we also need to include the prophylatic mastectomy group in the box and whisker plots.
#11May:Remove RNA140128LC_155, not marked as PM 
data1=read.table("mRNA_lowcov_removed_quantile_normalized_batch_corrected.txt",header=T,sep="\t")
annot=read.xlsx("annotation_13april2017.xlsx",1)
tumors=subset(annot, annot$tissue=="Tumor" | annot$tissue=="ProphylacticMasectomy"  )
tumors=subset(tumors,tumors$confidance!=1 )
tumors=subset(tumors,tumors$sample!="RNA140128LC_155" )

#TAM markers
markers=read.xlsx("tam_markers_revised.xlsx",1)
tamdata=data1[rownames(data1) %in% markers$tam,colnames(data1) %in% tumors$sample]

annotation_col = data.frame(tumors[,c(2,6,7)])
rownames(annotation_col)=colnames(tamdata)

subtype_order=order(tumors$subtype)
tamdata_order=tamdata[subtype_order]

#median normalize tam data
middle=apply(tamdata,1,median)
tamdata_n=tamdata-middle
tamdata_t=t(tamdata_n)
tamdata_t=as.matrix(tamdata_t)
tamscore=apply(tamdata_t,1,median)

tamscore=as.data.frame(tamscore)
tumors$tamscore=tamscore$tamscore

#the level order is changed 
tumors$subtype=gsub("ProphylacticMasectomy","PM",tumors$subtype)
tumors$subtype=factor(tumors$subtype,levels=c("Basal","Her2","LumA","LumB","TumorNormalLike","PM"))

# Let's get the sizes right, based on sample size 
levelProportions=summary(tumors$subtype)/50

pdf("tam_subtype_11may2017.pdf")
boxplot(tumors$tamscore~tumors$subtype,col=c("red","pink","darkblue","lightblue","green4","green3"),width=levelProportions,ylab="TAM signature score",xlab="BRCA subtypes")

# Add data points
mylevels=levels(tumors$subtype)
for(i in 1:length(mylevels))
{
  thislevel=mylevels[i]
  thisvalues=tumors[tumors$subtype==thislevel, "tamscore"]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter=jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.2)) 
  
  #  add  text
  TopOfWhisker=min(max(thisvalues), median(thisvalues)+IQR(thisvalues)*1.5)
  text(i+levelProportions[i]/2, TopOfWhisker, labels=paste("N=", length(thisvalues), sep=""), cex=.6, font=2, pos=4)
}
dev.off()

#ANOVA
summary(aov(tumors$tamscore~tumors$subtype))
aov1=aov(tumors$tamscore~tumors$subtype)
posthoc = TukeyHSD(x=aov1)
table(tumors$subtype)

