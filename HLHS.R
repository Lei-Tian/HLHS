# Combine all (LVn+RVn+LVp+RVp) cells of undeleoped left heart patient's fetal heart
library(Seurat)
library(dplyr)
setwd("D:/SharedFolder/project/Mingxia/Seurat")
ind="LVp"
load(paste0(ind,"/",ind,".Robj"))
LVp=myData
ind="LVn"
load(paste0(ind,"/",ind,".Robj"))
all=MergeSeurat(object1=LVp,object2=myData,add.cell.id1="LVp",add.cell.id2="LVn",do.normalize = FALSE)

for(ind in c("RVp","RVn")){
  load(paste0(ind,"/",ind,".Robj"))
  all=MergeSeurat(object1=all,object2=myData,add.cell.id2=ind,do.normalize = FALSE)
}

all <- FilterCells(object = all, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.2))
all <- NormalizeData(object = all, normalization.method = "LogNormalize", scale.factor = 10000)
pdf("HLHS/variable.genes.pdf")
all <- FindVariableGenes(object = all, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = all@var.genes)
all <- ScaleData(object = all, vars.to.regress = c("nUMI", "percent.mito"))
all <- RunPCA(object = all, pc.genes = all@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5, pcs.compute = 30)
pdf("HLHS/PCA.pdf")
VizPCA(object = all, pcs.use = 1:4)
VizPCA(object = all, pcs.use = 5:8)
PCAPlot(object = all, dim.1 = 1, dim.2 = 2)
dev.off()
all <- ProjectPCA(object = all, do.print = FALSE)
pdf("HLHS/PCA.heatmap.pdf")
PCHeatmap(object = all, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()
all <- JackStraw(object = all, num.replicate = 100, display.progress = FALSE, num.pc = 30)
pdf("HLHS/determine.NO.PCs.pdf")
JackStrawPlot(object = all, PCs = 1:30)
PCElbowPlot(object = all,num.pc =30)
dev.off()

all <- FindClusters(object = all, reduction.type = "pca", dims.use = 1:20, resolution = 0.8,force.recalc=T, print.output = 0, save.SNN = TRUE)
all <- RunTSNE(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)
all <- RunUMAP(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)

pdf("HLHS/tSNE.pdf")
DimPlot(all, reduction.use = "tsne", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "tsne", pt.size =0.6, group.by = "orig.ident")
dev.off()

pdf("HLHS/uMAP.pdf")
DimPlot(all, reduction.use = "umap", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "umap", pt.size =0.6, group.by = "orig.ident")
dev.off()

for(ct in c("CM","FB","SMC")){
  GeneList=read.table(paste0(ct,".markers.txt"))
  GeneList=toupper(GeneList[,1])
  GeneList=GeneList[GeneList %in% row.names(all@data)]
  for(gene in GeneList){
    print(gene)
    p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5,do.return =T)
    ggsave(p[[1]],filename = paste0("HLHS/uMAP/",ct,"/",gene,".pdf"))
  }
}

GeneList=c("CDH5","LYVE1","UPK3B","CD14","AHSP","S100B","CHN4")
for(gene in GeneList){
  print(gene)
  p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5,do.return =T)
  ggsave(p[[1]],filename = paste0("HLHS/uMAP/",gene,".pdf"))
}

for(gene in c("NPR3", "HAPLN1", "CDH11", "MGLL", "APLN", "APLNR")){
  print(gene)
  p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5, do.return =TRUE)
  ggsave(p[[1]],filename = paste0("HLHS/uMAP/20190911/",gene,".pdf"))
}

combined.markers <- FindAllMarkers(object = all,  min.pct = 0.25, thresh.use = log(2), only.pos = T)
combined.markers=combined.markers[which(combined.markers$p_val_adj<0.05),]
write.table(combined.markers,"HLHS/Clusters.DEGs.txt",quote=F,sep="\t")
top10 <- combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("HLHS/Clusters.DEG.heatmap.pdf",height=12,width=9)
DoHeatmap(object = all, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20)
dev.off()

library(RColorBrewer)
RdBu = rev(brewer.pal(11, name="RdBu"))
load("HLHS/all.Robj")
genes=c("TFE3","EDNRA","ZNF292","FOXM1","ZMYND19","PCBP3","TCF12","ARID1B","NOVA1","CACNA1",
        "PKD1","RBFOX2","ST5","TSC1","USP8","HERC4","KMT2D","ETS1","CHD7","CTR9","GLA",
        "FMNL1","PHRF1","SIPA1L1","HIRA")

genes=genes[genes %in% row.names(all@data)]
aveExp=c()
clusters=list("Endothelium"=c(0,4,5,8,15),"Endocarduim"=c(2,3,13,17),"Valve EC"=c(18),
              "Lymphatic EC"=c(19),"Immune Cell"=c(12,14),"HSC"=c(6,16),"Neuron"=c(21),
              "Fibrolast"=c(1,7,10,11),"SMC"=c(9),"CM"=c(20),"Epicarduim"=c(22))

for(ct in names(clusters)){
  print(ct)
  cells=names(all@ident)[which(all@ident %in% clusters[[ct]])]
  expression=all@data[genes,cells]
  meanExp=apply(expression,1,mean)
  aveExp=cbind(aveExp,meanExp)
}
colnames(aveExp)=names(clusters)
#aveExp=aveExp[-1,]
write.table(aveExp,"HLHS/HLHS.heatmap.merged.txt",sep="\t",quote=F)
exp.df=aveExp#[genes,c(6,9,12,13,1,11,22,17,14,20,0,2,3,4,7,15,19,8,10,5,16,18,21)]
library(dendsort)
library(gplots)
pdf("HLHS/HLHS.heatmap.merged.pdf",width=10,height=8)
heatmap.2(exp.df,Rowv=F,Colv=F,
          col=RdBu,srtCol=45,adjCol = c(1,1),
          dendrogram ="none",
          scale="row", labRow=rownames(exp.df), labCol=colnames(exp.df),
          cexRow=1.2, cexCol=1.2, margins=c(7, 7),
          symm = F, key = T, keysize =1, trace="none", density.info="none",lhei = c(1,4))
dev.off()
