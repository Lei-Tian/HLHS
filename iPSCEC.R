# Combine HLHS and normal iPSC-ECs
library(Seurat)
library(dplyr)
setwd("D:/SharedFolder/project/Mingxia/Seurat")

load("iPSCEC/all.Robj")

#all <- FindClusters(object = all, reduction.type = "pca", dims.use = 1:20, resolution = 0.8,force.recalc=T, print.output = 0, save.SNN = TRUE)
#all <- RunTSNE(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)
all <- RunUMAP(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F,max.dim=3)


cbbPalette <- c("#3399CC", "#339999", "#339966", "#993333", "#993366",  "#999933", "#CC6666")

pdf("iPSCEC/tSNE.pdf")
DimPlot(all, reduction.use = "tsne", do.label=T, pt.size = 1.2, label.size = 10)+
      theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
        legend.text = element_text( size = 20, face = "bold"))
DimPlot(all, reduction.use = "tsne", pt.size =1.2, group.by = "orig.ident")+
      theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
        legend.text = element_text( size = 20, face = "bold"))
dev.off()

pdf("iPSCEC/uMAP.pdf")
DimPlot(all, reduction.use = "umap", do.label=F, pt.size = 1.2, label.size = 10)+
  theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
        legend.text = element_text( size = 20, face = "bold"))
DimPlot(all, reduction.use = "umap", do.label=T, pt.size = 1.2, label.size = 10)+
      theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
        legend.text = element_text( size = 20, face = "bold"))
DimPlot(all, reduction.use = "umap", pt.size =1.2, group.by = "orig.ident")+
      theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
        legend.text = element_text( size = 20, face = "bold"))
dev.off()



markers=read.table("tSNE_plot_genes.txt")
markers=as.vector(markers[,1])
markers=c(markers,c("FN1","NTS","APOE","COLEC11","IGFBP5","C7","NPR3","APCDD1","CFH","ENG","CGNL1","SCD","CCDC80","IGFBP3","TMEM100","LDB2","TMSB4X","SERPINE2","ADGRG6"))
markers=markers[markers %in% row.names(all@data)]

markers=c("FN1","CGNL1","ADGRG6","TFPI","SAT1","NPR3","HAPLN1","PLVAP","CDH11","MGLL","APLN","TGFB1","TGFB2","SNAI2","FN1","TAGLN","NOTCH1","DLL4","JAG1","HEY1")
for(gene in markers){
  print(gene)
  p=FeatureHeatmap(object = all, features.plot = c(gene), group.by="orig.ident", reduction.use = "umap",cols.use=c("lightgrey", "blue"),pt.size=1.2,plot.horiz=T)+
    theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
          legend.text = element_text( size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("iPSCEC/uMAP/",gene,".pdf"),height=7.5,width=6)
}

all@meta.data$tmp <- factor(x = all@ident, levels = c(2,0,5,3,6,4,1))
cbbPalette <- c("#339966", "#3399CC", "#999933", "#993333", "#CC6666", "#993366", "#339999")
for(gene in markers){
  print(gene)
  p=VlnPlot(object = all, features.plot = c(gene),group.by="tmp",point.size.use = 0,cols.use =cbbPalette)+
    geom_boxplot(width=.05,outlier.alpha=0,fill="black",colour="gray",coef=0.1) + 
    theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
          axis.text.x  = element_text( vjust=0.5, size=34,face="bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  #+VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  #p2=VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  ggsave(p,filename = paste0("iPSCEC/Final Violin/",gene,".pdf"),height=6,width=12)
}

save(all,file="iPSCEC/all.Robj")

c026=SubsetData(all, ident.use = c(0,2,6))
c1345=SubsetData(all, ident.use = c(1,3,4,5)) 
c026=SetAllIdent(c026,id="orig.ident")
c1345=SetAllIdent(c1345,id="orig.ident")

c026.markers <- FindMarkers(object = c026, ident.1="M146", ident.2 = "BA060")
write.table(c026.markers,"iPSCEC/c026.M146vsBA060.txt",quote=F,sep="\t")
c1345.markers <- FindMarkers(object = c1345, ident.1="M146", ident.2 = "BA060")
write.table(c1345.markers,"iPSCEC/c1345.M146vsBA060.txt",quote=F,sep="\t")


genes=read.table("iPSCEC/genes_mutation_expression.txt")
genes=as.vector(genes[,1])
pdf("iPSCEC/mutated.genes.heatmap.pdf",height=12,width=9)
DoHeatmap(object = all, genes.use = genes, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20)
DoHeatmap(object = all, genes.use = genes, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20, group.by = "orig.ident")
dev.off()
all=SetAllIdent(all,id="orig.ident")
combined.markers <- FindMarkers(object = all, genes.use =genes, ident.1="M146", ident.2 = "BA060", min.pct = 0, logfc.threshold = 0)
write.table(combined.markers,"iPSCEC/mutated.genes.statistics.txt",quote=F,sep="\t")


kkk=sample(which(tmp$group=="others"),25)
tmp[kkk,"group"]="random"
bbb=tmp[which(tmp$group %in% c("Mutation","random")),]
t.test(bbb[which(bbb$group=="Mutation"),"avg_logFC"],bbb[which(bbb$group=="Random"),"avg_logFC"])
pdf("random.pdf",height=7,width=8)
par(mar=c(5,8,5,8))
ggplot(bbb,aes(x=group,y=avg_logFC,colour=group))+geom_boxplot(alpha=0.8)+geom_point(size=3)+
  ylab(bquote('log'['2']*'FC(HLHS/Control)'))+
  theme(axis.text.x = element_text(size=28),axis.text.y = element_text(size=28),axis.title.y = element_text(size=32,face="bold"),legend.position="none")+
  scale_x_discrete(name="")+ ggtitle("iPSC-EC: HLHS v.s. Control")+scale_fill_manual(values=c("red", "blue"))+scale_colour_manual(values=c("red", "blue"))+
  theme(plot.title = element_text(face="bold",size=36,colour="#990000",hjust=1.1))
dev.off()

genes=c("HAPLN1","PLVAP","FN1","LIMCH1","JUNB","PXDN","CD24","PALMD","SPTBN1","TGM2","NREP","LMO2","LDB2","MAP4K4","TM4SF1","CD9")
gene_order=c("FN1","PXDN","SPTBN1","TGM2","LDB2","MAP4K4","LIMCH1","JUNB","HAPLN1","PLVAP","CD24","PALMD","LMO2","NREP","TM4SF1","CD9")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(0,2)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster==6),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endocarduim.FC.barplot.pdf",height=9,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="forestgreen",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(gene_order))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-2.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

genes=c("CALM2","SRGN","S100A10","HMGB1","CD59","ENO1","CAV1","MT2A","TAGLN","HIST1H4C")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(3,5)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster %in% c(1,4)),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endothelium.FC.barplot.pdf",height=6.5,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="orange",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=degs1$gene_name)+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-2,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

load("iPSCEC/all.Robj")
genes=c("UBE2S","SMC4","MT2A","CAV1","CD59","HMGB1","SRGN","CALM2")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(3,5)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster %in% c(1,4)),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endothelium.FC.barplot.RV.pdf",height=6.5,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="orange",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(genes))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

genes=c("PALMD","HAPLN1","JUNB","CD24","ITM2C","TM4SF1","CD9","HLA-C")
gene_order=c("PALMD","HAPLN1","JUNB","CD24","ITM2C","TM4SF1","CD9","HLA-C")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(0,2)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster==6),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endocarduim.FC.barplot.RV.pdf",height=9,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="forestgreen",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(gene_order))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

load("iPSCEC/all.Robj")
newID=rep("c0235",nrow(all@meta.data))
newID[which(all@ident %in% c(1,4,6))]="c146"
all@meta.data$newID=newID
all=SetAllIdent(all,id="newID")
c0235.vs.c146 <- FindMarkers(object = all, ident.1="c0235", ident.2 = "c146")
tmp <- FindMarkers(object = all,genes.use = c("ETS1","NOTCH1","DLL4", "JAG1", "HEY1"), ident.1="c0235", ident.2 = "c146", min.pct=0, logfc.threshold = 0)
c0235.vs.c146=rbind(c0235.vs.c146,tmp)
write.table(c0235.vs.c146,"iPSCEC/c0235.vs.c146.txt",quote=F,sep="\t")
