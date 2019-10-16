# Combine all (LVn+RVn+LVp+RVp) cells of normal fetal heart

library(Seurat)
library(dplyr)
setwd("D:/SharedFolder/project/Mingxia/Seurat")

ind="LVp"
load(paste0("E:/Mingxia/v2/Seurat/d83_",ind,"/","d83_",ind,".Robj"))
LVp=myData
ind="LVn"
load(paste0("E:/Mingxia/v2/Seurat/d83_",ind,"/","d83_",ind,".Robj"))
all=MergeSeurat(object1=LVp,object2=myData,add.cell.id1="LVp",add.cell.id2="LVn",do.normalize = FALSE)

for(ind in c("RVp","RVn")){
  load(paste0("E:/Mingxia/v2/Seurat/d83_",ind,"/","d83_",ind,".Robj"))
  all=MergeSeurat(object1=all,object2=myData,add.cell.id2=ind,do.normalize = FALSE)
}

all <- FilterCells(object = all, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.2))

all <- NormalizeData(object = all, normalization.method = "LogNormalize", scale.factor = 10000)
pdf("d83/variable.genes.pdf")
all <- FindVariableGenes(object = all, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = all@var.genes)
all <- ScaleData(object = all, vars.to.regress = c("nUMI", "percent.mito"))
all <- RunPCA(object = all, pc.genes = all@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5, pcs.compute = 30)
pdf("d83/PCA.pdf")
VizPCA(object = all, pcs.use = 1:4)
VizPCA(object = all, pcs.use = 5:8)
PCAPlot(object = all, dim.1 = 1, dim.2 = 2)
dev.off()
all <- ProjectPCA(object = all, do.print = FALSE)
pdf("d83/PCA.heatmap.pdf")
PCHeatmap(object = all, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()
all <- JackStraw(object = all, num.replicate = 100, display.progress = FALSE, num.pc = 30)
pdf("d83/determine.NO.PCs.pdf")
JackStrawPlot(object = all, PCs = 1:30)
PCElbowPlot(object = all,num.pc =30)
dev.off()

all <- FindClusters(object = all, reduction.type = "pca", dims.use = 1:20, resolution = 0.8,force.recalc=T, print.output = 0, save.SNN = TRUE)
all <- RunTSNE(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)
all <- RunUMAP(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)

pdf("d83/tSNE.pdf")
DimPlot(all, reduction.use = "tsne", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "tsne", pt.size =0.6, group.by = "orig.ident")
dev.off()

pdf("d83/uMAP.pdf")
DimPlot(all, reduction.use = "umap", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "umap", pt.size =0.6, group.by = "orig.ident")
dev.off()

for(ct in c("CM","FB","SMC")){
  GeneList=read.table(paste0(ct,".markers.txt"))
  GeneList=toupper(GeneList[,1])
  GeneList=GeneList[GeneList %in% row.names(all@data)]
  for(gene in GeneList){
    print(gene)
    p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5, do.return =TRUE)
    ggsave(p[[1]],filename = paste0("d83/uMAP/",ct,"/",gene,".pdf"))
  }
}

for(gene in c("CDH5","UPK3B","WT1","KRT8","KRT19","PROX1","LYVE1","CCL21","AHSP","SLC4A1","NFE2","GATA1","CD14",
              "CD68","KLRB1","NKG7","CD37","ALOX5AP","FCER1G","FST","PLP1","MPZ","HCN4","TBX3","TBX5")){
  print(gene)
  p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5, do.return =TRUE)
  ggsave(p[[1]],filename = paste0("d83/uMAP/20190404/",gene,".pdf"))
}

for(gene in c("NPR3", "HAPLN1", "CDH11", "MGLL", "APLN", "APLNR")){
  print(gene)
  p=FeaturePlot(object = all, features = c(gene), reduction = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5, do.return =TRUE)
  ggsave(p[[1]],filename = paste0("d83/uMAP/20190911/",gene,".pdf"))
}
markers=read.table("tSNE_plot_genes.txt")
markers=as.vector(markers[,1])
markers=markers[markers %in% row.names(all@data)]

for(gene in markers){
  print(gene)
  p=FeatureHeatmap(object = all, features.plot = c(gene), group.by="orig.ident", reduction.use = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5,plot.horiz=T)
  ggsave(p,filename = paste0("d83/uMAP/",gene,".pdf"),height=18,width=6)
}


combined.markers <- FindAllMarkers(object = all,  min.pct = 0.25, thresh.use = log(2), only.pos = T)
combined.markers=combined.markers[which(combined.markers$p_val_adj<0.05),]
write.table(combined.markers,"d83/Clusters.DEGs.txt",quote=F,sep="\t")
top10 <- combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("d83/Clusters.DEG.heatmap.pdf",height=12,width=9)
DoHeatmap(object = all, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20)
dev.off()

load("d83/all.Robj")
all <- StashIdent(object = all, save.name = "cluster")
all@meta.data$orig="RV"
all@meta.data$orig[which(all@meta.data$orig.ident %in% c("d83_LVp","d83_LVn"))]="LV"
all@meta.data$orig.cluster=paste0(all@meta.data$orig,"_",all@ident)
all <- SetAllIdent(all,id="orig.cluster")
results=c()
for(cluster in names(table(all@meta.data$cluster))){
  print(cluster)
  tmp <- FindMarkers(object = all, ident.1 = paste0("LV_",cluster), ident.2 = paste0("RV_",cluster))
  gene_name=row.names(tmp)
  tmp=cbind(tmp,gene_name,cluster)
  results=rbind(results,tmp)
}
write.table(results,"d83/LV_vs_RV.clusters.DEGs.txt",quote=F,sep="\t")


save(all,file="d83/all.Robj")

library(RColorBrewer)
RdBu = rev(brewer.pal(11, name="RdBu"))
load("d83/all.Robj")
genes=c("TFE3","EDNRA","ZNF292","FOXM1","ZMYND19","PCBP3","TCF12","ARID1B","NOVA1","CACNA1",
        "PKD1","RBFOX2","ST5","TSC1","USP8","HERC4","KMT2D","ETS1","CHD7","CTR9","GLA",
        "FMNL1","PHRF1","SIPA1L1","HIRA")

genes=genes[genes %in% row.names(all@data)]
aveExp=c()
for(c in c(6,9,12,13,1,11,22,17,14,20,0,2,3,4,7,15,19,8,10,5,16,18,21)){
  print(c)
  cells=names(all@ident)[which(all@ident==c)]
  expression=all@data[genes,cells]
  meanExp=apply(expression,1,mean)
  aveExp=cbind(aveExp,c(c,meanExp))
}
colnames(aveExp)=aveExp[1,]
aveExp=aveExp[-1,]
write.table(aveExp,"d83/d83.heatmap.txt",sep="\t",quote=F)
exp.df=aveExp#[genes,c(6,9,12,13,1,11,22,17,14,20,0,2,3,4,7,15,19,8,10,5,16,18,21)]
#Endocarduim, Endothelium, HSC, Lym EC, Immune Cell, Nervous system, Fibo, SMC, CM, Epicarduim
library(dendsort)
exp.df=t(apply(exp.df,2,function(x){return(x/max(x))}))
pdf("d83/test.pdf",width=10,height=6)
heatmap.2(exp.df,Rowv=F,Colv=F,
          col=RdBu,srtCol=0,adjCol = c(0.5,1),
          dendrogram ="none",
          scale="row", labRow=rownames(exp.df), labCol=colnames(exp.df),
          cexRow=1.2, cexCol=1.5, margins=c(5, 7),
          symm = F, key = T, keysize =1, trace="none", density.info="none",lhei = c(1,4))
dev.off()

clusters=list("Endothelium"=c(6,9,12,13),"Endocarduim"=c(1,11),"HSC"=c(22),
              "Lymphatic EC"=c(17),"Immune Cell"=c(14),"Neuron"=c(20),
              "Fibrolast"=c(0,2,3,4,7,15,19),"SMC"=c(8,10),"CM"=c(5,16),
              "Conduction System"=c(18),"Epicarduim"=c(21))
aveExp=read.table("d83/d83.heatmap.txt",sep="\t",check.names = F)
#exp.df=as.matrix(aveExp)
maxExp=c()
for(ct in names(clusters)){
  print(ct)
  if(length(clusters[[ct]])>1){
    expression=apply(aveExp[,as.character(clusters[[ct]])],1,max)
  }else{
    expression=aveExp[,as.character(clusters[[ct]])]
  }
  maxExp=cbind(maxExp,expression)
}
colnames(maxExp)=names(clusters)
row.names(maxExp)=row.names(exp.df)
pdf("d83/d83.heatmap.max.pdf",width=10,height=8)
heatmap.2(maxExp,Rowv=F,Colv=F,
          col=RdBu,srtCol=45,adjCol = c(1,1),
          dendrogram ="none",
          scale="row", labRow=rownames(maxExp), labCol=colnames(maxExp),
          cexRow=1.2, cexCol=1.2, margins=c(7, 7),
          symm = F, key = T, keysize =1, trace="none", density.info="none",lhei = c(1,4))
dev.off()

load("d83/all.Robj")
genes=c("TFE3","EDNRA","ZNF292","FOXM1","ZMYND19","PCBP3","TCF12","ARID1B","NOVA1","CACNA1",
        "PKD1","RBFOX2","ST5","TSC1","USP8","HERC4","KMT2D","ETS1","CHD7","CTR9","GLA",
        "FMNL1","PHRF1","SIPA1L1","HIRA")

genes=genes[genes %in% row.names(all@data)]
aveExp=c()
clusters=list("Endothelium"=c(6,9,12,13),"Endocarduim"=c(1,11),"HSC"=c(22),
              "Lymphatic EC"=c(17),"Immune Cell"=c(14),"Neuron"=c(20),
              "Fibrolast"=c(0,2,3,4,7,15,19),"SMC"=c(8,10),"CM"=c(5,16),
              "Conduction System"=c(18),"Epicarduim"=c(21))
for(ct in names(clusters)){
  print(ct)
  cells=names(all@ident)[which(all@ident %in% clusters[[ct]])]
  expression=all@data[genes,cells]
  meanExp=apply(expression,1,mean)
  aveExp=cbind(aveExp,meanExp)
}
colnames(aveExp)=names(clusters)
write.table(aveExp,"d83/d83.heatmap.merged.txt",sep="\t",quote=F)
aveExp=read.table("d83/d83.heatmap.merged.txt",sep="\t")
exp.df=as.matrix(aveExp)
library(dendsort)
library(gplots)
RdBu = rev(brewer.pal(11, name="RdBu"))
pdf("d83/d83.heatmap.merged.non.pdf",width=10,height=8)
heatmap.2(exp.df,Rowv=F,Colv=F,
          col=RdBu,srtCol=45,adjCol = c(1,1),
          dendrogram ="none",
          scale="none", labRow=rownames(exp.df), labCol=colnames(exp.df),
          cexRow=1.2, cexCol=1.2, margins=c(7, 7),
          symm = F, key = T, keysize =1, trace="none", density.info="none",lhei = c(1,4))
dev.off()
