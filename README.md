# 20250208
#install.packages("pheatmap")
#install.packages("stringr")
#install.packages("ggplot2")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("devtools")
#devtools::install_github("BioSenior/ggVolcano")

library(limma)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggVolcano)

setwd("F:/pSS_HT")
data=read.table("GSE.txt", header=T, sep="\t", check.names=F,row.names = 1)
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#data=data[rowMeans(data)>1,]
boxplot(data.frame(data),col="#4DBBD5")
data=normalizeBetweenArrays(data)
boxplot(data.frame(data),col="#4DBBD5")
dev.off()
write.table(data.frame(ID=rownames(data),data),file="normalize.txt", sep="\t", quote=F, row.names = F)
Control=read.table("Control.txt", header=F, sep="\t", check.names=F)
Case=read.table("Caset.txt", header=F, sep="\t", check.names=F)
conNum=length(rownames(Control))
treatNum=length(rownames(Case))
Type=c(rep(1,conNum), rep(2,treatNum))
data1 = data[,Control[,1]]
data2 = data[,Case[,1]]
data = cbind(data1,data2)
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  #logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  logFC=treatGeneMeans-conGeneMeans
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)
write.table(outTab,file="all.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

logFCfilter=1.5
fdrFilter=0.05
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & 
                   as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

geneNum=50     
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=6.5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("#4DBBD5",5), "white", rep("#E64B35",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

pdf(file="vol.pdf", width=5, height=5)
xMax=6
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1.2)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#E64B35",cex=1.5)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#4DBBD5",cex=1.5)
abline(v=0,lty=2,lwd=3)
dev.off()

library(VennDiagram)
venn.diagram(list("GSE1" = List[[1]],
                  "GSE2" = List[[2]]) ,
             height=5000,
             width=5200,
             resolution=500,
             imagetype="tiff",
             filename="VennPlot2.tiff",
             col="transparent",
             fill = c("darkblue", "orange"),alpha = 0.50,
             label.col = c("orange", "white", "darkblue", "white"),
             cat.col = c("darkblue", "orange"),
             cex=1.5,
             cat.cex=1.4)
             
rt = read.table(file = 'pathway.txt',sep = '\t',header = T,quote = '')
keggSig = rt[rt$fdr < 0.05,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))

library(ggplot2)
ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
    color=expression(-log[10](fdr)),
    size="Gene number",
    x="Fold enrichment"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
ggsave('plot.pdf',width = 7,height = 6)

library(IOBR)
data=read.table("GSE.txt", header=T, sep="\t", check.names=F,row.names = 1)
#data = log2(data +1)
#data=data[rowMeans(data)>0.5,]
# method:'mcpcounter', 'epic', 'xcell', 'cibersort', 
im_mcpcounter <- deconvo_tme(eset = data,method = "mcpcounter")
im_epic <- deconvo_tme(eset = data,method = "epic",arrays = F)
im_xcell <- deconvo_tme(eset = data,method = "xcell",arrays = F)
im_cibersort <- deconvo_tme(eset = data,method = "cibersort",arrays = F,perm = 1000))
im_genesets = signature_collection
im_ssgsea <- calculate_sig_score(eset=data,signature=signature_collection,method= "ssgsea")
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
saveRDS(tme_combine,"tme_combine.rds")
tme_combine = readRDS('tme_combine.rds')
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
group=gsub("1", "Control", group)
group=gsub("0", "Case", group)
Type = cbind(tme_combine[,1],group)
signature_group = sig_group
names(signature_group)[1]
names(signature_group)[1] = "Case_Signature"
names(signature_group)[1]
iobr_cor_plot(pdata_group = Type,id1 = "ID",
              feature_data = tme_combine,id2 = "ID",
              group = "group",is_target_continuous  = F,
              category = "signature",character_limit=60,
              #signature_group = sig_group[c(1:10)],
              palette_box = "jco",ProjectID = "GSE")
