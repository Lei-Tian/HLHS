# The codes are combine the right ventricle CDH5+ cells from control and 
# Undeveloped Left Heart patient. And perform comparison
library(Seurat)
library(dplyr)
setwd("D:/SharedFolder/project/Mingxia/Seurat")
ind="RVp"
load(paste0("E:/Mingxia/v2/Seurat/d83_RVp/d83_RVp.Robj"))
d83=myData
load("RVp/RVp.Robj")
all=MergeSeurat(object1=d83,object2=myData,add.cell.id1="d83",add.cell.id2="HLHS",do.normalize = FALSE)

all <- FilterCells(object = all, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.2))

all <- NormalizeData(object = all, normalization.method = "LogNormalize", scale.factor = 10000)
pdf("d83_HLHS_RVp/variable.genes.pdf")
all <- FindVariableGenes(object = all, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = all@var.genes)
all <- ScaleData(object = all, vars.to.regress = c("nUMI", "percent.mito"))
all <- RunPCA(object = all, pc.genes = all@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5, pcs.compute = 30)
pdf("d83_HLHS_RVp/PCA.pdf")
VizPCA(object = all, pcs.use = 1:4)
VizPCA(object = all, pcs.use = 5:8)
PCAPlot(object = all, dim.1 = 1, dim.2 = 2)
dev.off()
all <- ProjectPCA(object = all, do.print = FALSE)
pdf("d83_HLHS_RVp/PCA.heatmap.pdf")
PCHeatmap(object = all, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()
all <- JackStraw(object = all, num.replicate = 100, display.progress = FALSE, num.pc = 30)
pdf("d83_HLHS_RVp/determine.NO.PCs.pdf")
JackStrawPlot(object = all, PCs = 1:30)
PCElbowPlot(object = all,num.pc =30)
dev.off()

all <- FindClusters(object = all, reduction.type = "pca", dims.use = 1:20, resolution = 0.8,force.recalc=T, print.output = 0, save.SNN = TRUE)
all <- RunTSNE(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)
all <- RunUMAP(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)

pdf("d83_HLHS_RVp/tSNE.pdf")
DimPlot(all, reduction.use = "tsne", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "tsne", pt.size =0.6, group.by = "orig.ident")
dev.off()

pdf("d83_HLHS_RVp/uMAP.pdf")
DimPlot(all, reduction.use = "umap", pt.size =0.6,do.label=T)
DimPlot(all, reduction.use = "umap", pt.size =0.6, group.by = "orig.ident")
dev.off()


combined.markers <- FindAllMarkers(object = all,  min.pct = 0.25, thresh.use = log(2), only.pos = T)
combined.markers=combined.markers[which(combined.markers$p_val_adj<0.05),]
write.table(combined.markers,"d83_HLHS_RVp/Clusters.DEGs.txt",quote=F,sep="\t")
top10 <- combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("d83_HLHS_RVp/Clusters.DEG.heatmap.pdf",height=12,width=9)
DoHeatmap(object = all, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20)
dev.off()

markers=read.table("tSNE_plot_genes.txt")
markers=as.vector(markers[,1])
markers=markers[markers %in% row.names(all@data)]

for(gene in markers){
  print(gene)
  p=FeatureHeatmap(object = all, features.plot = c(gene), group.by="orig.ident", reduction.use = "umap",cols.use=c("lightgrey", "blue"),pt.size=0.5,plot.horiz=T)
  ggsave(p,filename = paste0("d83_HLHS_RVp/uMAP/",gene,".pdf"),height=7.5,width=6)
  p=VlnPlot(object = all, features.plot = c(gene),point.size.use = 0)+
    geom_boxplot(width=.05,outlier.alpha=0,fill="black",colour="gray",coef=0.1) + 
    theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
          axis.text.x  = element_text( vjust=0.5, size=24,face="bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("d83_HLHS_RVp/violin/",gene,".pdf"),height=6,width=12)
}


all <- StashIdent(object = all, save.name = "cluster")
all@meta.data$orig.cluster=paste0(all@meta.data$orig.ident,"_",all@ident)
all <- SetAllIdent(all,id="orig.cluster")
results=c()
for(cluster in names(table(all@meta.data$cluster))){
  print(cluster)
  tmp <- FindMarkers(object = all, ident.1 = paste0("d83_RVp_",cluster), ident.2 = paste0("RVp_",cluster))
  gene_name=row.names(tmp)
  tmp=cbind(tmp,gene_name,cluster)
  results=rbind(results,tmp)
}
write.table(results,"d83_HLHS_RVp/d83_vs_HLHS.clusters.DEGs.txt",quote=F,sep="\t")

write.table(table(all@ident),"d83_HLHS_RVp/cell.numers.txt",quote=F,sep="\t")

save(all,file="d83_HLHS_RVp/all.Robj")

genes=c("PALMD","HAPLN1","JUNB","CD24","ITM2C","TM4SF1","CD9","HLA-C")
gene_order=genes
degs=read.table("d83_HLHS_RVp/d83_vs_HLHS.clusters.DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene_name %in% genes),]
degs=degs[which(degs$cluster %in% c(1,2,5,6,9,10)),]
degs$avg_logFC=-degs$avg_logFC
degs$absFC=abs(degs$avg_logFC)
degs1=degs %>% group_by(gene_name) %>% top_n(1,absFC)
degs1=degs1[order(degs1$avg_logFC,decreasing = TRUE),]
pdf("d83_HLHS_RVp/d83_HLHS_RVp.Endocarduim.FC.barplot.pdf",height=9,width=6)
ggplot(degs1,aes(x=gene_name,y=avg_logFC))+geom_col(fill="forestgreen",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(gene_order))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

genes=c("UBE2S","SMC4","MT2A","CAV1","CD59","HMGB1","SRGN","CALM2")
degs=read.table("d83_HLHS_RVp/d83_vs_HLHS.clusters.DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene_name %in% genes),]
degs=degs[which(degs$cluster %in% c(0,3,4,7,8)),]
degs$avg_logFC=-degs$avg_logFC
degs$absFC=abs(degs$avg_logFC)
degs1=degs %>% group_by(gene_name) %>% top_n(1,absFC)
degs1=degs1[order(degs1$avg_logFC,decreasing = TRUE),]
pdf("d83_HLHS_RVp/d83_HLHS_RVp.Endothelium.FC.barplot.pdf",height=6.5,width=6)
ggplot(degs1,aes(x=gene_name,y=avg_logFC))+geom_col(fill="orange",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(genes))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()
