#! /usr/bin/env Rscript

library(DESeq2)
library(genefilter)
library(ggplot2)

source("/labs/joewu/tianlei/projects/lamc/RNAseq/my_PCA.R")

exLim <- function(values, ratio) {
  lim.min <- min(values)
  lim.max <- max(values)
  lim.mid <- (lim.max-lim.min)/2 + lim.min
  lim.min.new <- lim.mid - (lim.max-lim.min)*(1+ratio)/2
  lim.max.new <- lim.mid + (lim.max-lim.min)*(1+ratio)/2
  lims <- list(min=lim.min.new, max=lim.max.new)
  lims
}

FC.cutoff <- 2.0
args=commandArgs(T)
tag <- "HLHS"
#setwd(tag)
reads.cnt.tbl <- read.table("A.out.txt",
                            stringsAsFactors=FALSE,
                            header=TRUE, sep="\t")

rownames(reads.cnt.tbl) <- reads.cnt.tbl[ ,1]
reads.cnt.tbl <- reads.cnt.tbl[ , -1]
maxCounts=apply(reads.cnt.tbl,1,max)
reads.cnt.tbl=reads.cnt.tbl[which(maxCounts>=2),]

n=ncol(reads.cnt.tbl)
phenotypes=c(rep("CTR",3),rep("HLHS",3))
myCols=c(rep("blue",3),rep("red",3))
colData <- as.data.frame(cbind(phenotypes=phenotypes))
rownames(colData) <- names(reads.cnt.tbl)

print(colData)

cds <- DESeqDataSetFromMatrix(countData = reads.cnt.tbl,
                              colData = colData,
                              design = ~ phenotypes+1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds)
vsd.exp <- assay(vsd)
write.table(vsd.exp, file="vsd_exp.txt", sep="\t", quote=FALSE)
#
#cds <- DESeq(cds)
#print(resultsNames(cds))
#
set.seed(5)

pdf(file=paste("PCA_", tag, "_DESeq2.pdf", sep=""), width=8, height=8)

ntop=500
rv <- rowVars(vsd.exp)
select=order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
#select <- which(rv>0.1)#order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
vsd.exp=vsd.exp[select,]
pca <- prcomp(t(vsd.exp))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
write.table(pca$x,"pcs.txt",quote=F)
par(mar=c(5,5,2,2))
plot(pca$x[,"PC1"], pca$x[,"PC2"],cex=4, xlab = paste("PC1 (",round(100*percentVar[1],2),"%)",sep=""),ylab =paste("PC2 (",round(100*percentVar[2],2),"%)",sep=""),pch=16,cex.axis=2,cex.lab=2,col=myCols)
text(pca$x[,"PC1"], pca$x[,"PC2"],labels = colnames(vsd.exp))
#legend("bottomright",legend=c("K3-6","K1-2"),col=c("red","blue"),pch=16,cex=4)

plot(pca$x[,"PC1"], pca$x[,"PC3"],cex=4, xlab = paste("PC1 (",round(100*percentVar[2],2),"%)",sep=""),ylab =paste("PC3 (",round(100*percentVar[3],2),"%)",sep=""),pch=16,cex.axis=2,cex.lab=2,col=myCols)
text(pca$x[,"PC1"], pca$x[,"PC3"],labels = colnames(vsd.exp))
#legend("topright",legend=c("K3-6","K1-2"),col=c("red","blue"),pch=16,cex=4)

plot(pca$x[,"PC2"], pca$x[,"PC3"],cex=4, xlab = paste("PC2 (",round(100*percentVar[2],2),"%)",sep=""),ylab =paste("PC3 (",round(100*percentVar[3],2),"%)",sep=""),pch=16,cex.axis=2,cex.lab=2,col=myCols)
text(pca$x[,"PC2"], pca$x[,"PC3"],labels = colnames(vsd.exp))
#legend("topleft",legend=c("K5.6","K1-2"),col=c("red","blue"),pch=16,cex=4)

dev.off()

#print(colData(cds))
cds.1 <- DESeq(cds, test="LRT",full=~phenotypes+1,reduced=~1)
print(resultsNames(cds.1))
res.LRT.1 <- results(cds.1,name="phenotypes_HLHS_vs_CTR")
print(mcols(res.LRT.1)$description)
write.table(res.LRT.1, file="res_LRT.txt", sep="\t", quote=FALSE)


