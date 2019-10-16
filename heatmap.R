# library(gapmap)
library(dendsort)
library("gplots")
library(RColorBrewer)
#library(iterpc)
col.cut.k <- 2
row.cut.k <- 2
args=commandArgs(T)

RdBu = rev(brewer.pal(11, name="RdBu"))
RdYlBu = rev(brewer.pal(11, name="RdYlBu"))

exp.df <- read.table("DEGs.exp.txt", stringsAsFactors=F, header=T, row.names=1,check.names = F)

my.heatmap <- function(x, exp.df, row.cut.k, col.cut.k){
  print(x)
  cor.method <- x[1]
  hclust.method <- x[2]
  dendsort.method <- x[3]
  hc.r <- hclust(as.dist(1-cor(t(exp.df), method=cor.method)), method=hclust.method)
  gr.r <- cutree(hc.r, k=row.cut.k)
  hc.c <- hclust(as.dist(1-cor(exp.df, method=cor.method)), method=hclust.method)
  gr.c <- cutree(hc.c, k=col.cut.k)
  dend.r <- as.dendrogram(hc.r)
  dend.c <- as.dendrogram(hc.c)
  heatmap.2(as.matrix(exp.df),
          col=RdBu,
          dendrogram ="both",
           Rowv=rev(dendsort(dend.r, type = dendsort.method)), Colv=dendsort(dend.c, type = dendsort.method),
           scale="row", labRow=rownames(exp.df),
           labCol=colnames(exp.df),
           cexRow=2, 
           cexCol=2.5, margins=c(5, 15),
           RowSideColors=RdYlBu[gr.r], ColSideColors=c(rep("blue",3),rep("red",3)),
           symm = F, key = T, keysize =1, trace="none", density.info="none",lhei = c(1,7))#,
           #main=paste0("1-cor, ", cor.method, " ", hclust.method, " ", dendsort.method))
}

cor.methods <- c("spearman")
hclust.methods=c("ward.D2")
dendsort.methods <- c("average")
comb.vars <- expand.grid(cor.methods, hclust.methods, dendsort.methods)

pdf(file="DEGs.exp.pdf", width=6, height=10)
  res <- apply(comb.vars, 1, my.heatmap, exp.df=exp.df, row.cut.k=row.cut.k, col.cut.k=col.cut.k)
dev.off()

