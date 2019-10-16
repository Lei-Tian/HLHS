#! /usr/bin/env Rscript

# Usage: Rscript enrich_pathway_geneAnswer.R tag gene_list.yaml
# Example of yaml:
#     down: down_genes.txt
#     up: up_genes.txt
library(KEGG.db)
library(reactome.db)
library(GO.db)
library(yaml)
library(org.Hs.eg.db)
library(GeneAnswers)
library(ggplot2)
library(reshape2)
library(plyr)

e2s = toTable(org.Hs.egSYMBOL)
config <- yaml.load_file("enrich_geneAnswer.yaml")
args <- commandArgs(trailingOnly = TRUE)
tag <- "HLHS"
#setwd(tag)
species.anno <- "org.Hs.eg.db"
pvalueT <- 1

readGeneList <- function(gene.list.file.name, e2s){
  genes <- read.table(gene.list.file.name, header=FALSE, stringsAsFactors=FALSE)[,1]
  genes.entrez <- e2s$gene_id[e2s$symbol %in% genes]
}

output.genes <- function(gAn.ins, tag, categoryType){
  cat.funs <- c(GO="topGOGenes", DOLITE="topDOLITEGenes", KEGG="topPATHGenes",
    reactome.path="topREACTOME.PATHGenes")
  # cat.funs <- as.character(data.frame(cat=c("GO", "DOLITE", "KEGG", "reactome.path"),
  #   funs=c("topGOGenes", "topDOLITEGenes", "topPATHGenes", "topREACTOME.PATHGenes")))
  print(cat.funs)
  fun.name <- cat.funs[categoryType]
  FUN <- match.fun(fun.name)
  pathways.num <- length(gAn.ins@genesInCategory)
  genes.num=max(unlist(llply(gAn.ins@genesInCategory,length)))
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  write.table(gAn.ins@enrichmentInfo,paste0(tag,"_", sub.name,"_", categoryType, "_enrichmentInfo.txt"),quote=F,sep="\t")
  FUN(gAn.ins, orderby='pvalue', top=pathways.num,topGenes=genes.num,
    fileName=paste0(tag,"_", sub.name,"_", categoryType, "_genes.txt"), file=TRUE)
}

gAnAnalyze.lev3 <- function(categoryType, tag, genes.list, species.anno, pvalueT){
    gAn <- lapply(genes.list, geneAnswersBuilder,
                       species.anno, categoryType=categoryType, pvalueT=pvalueT,
                       level=3)
    print(categoryType)
    for(i in 1:length(gAn)){
      attr(gAn[[i]],"name") <- names(gAn)[i]
    }
    if(categoryType %in% c("GO", "DOLITE", "KEGG", "reactome.path")){
      try.out.genes <- llply(gAn, output.genes, tag=tag, categoryType=categoryType)
    }

    gAn.cluster <- getConceptTable(gAn, items='geneNum',topCat=10000000000)
    write.table(-log10(gAn.cluster$IndexTable),
      file=paste0(tag, "_enrich_", categoryType, "_lev3.txt"),
      sep="\t", quote=FALSE, row.names=TRUE)

    g.N <- gAn.cluster$CategoriesTable
    g.N$terms <- row.names(g.N)

    ## keep the order of the terms
    g.N$terms <- factor(g.N$terms, levels=rev(g.N$terms))

    g.N <- g.N[!g.N$terms=="Genes / Group", ]
    g.N.long <- melt(g.N)
    terms.p <- as.data.frame(gAn.cluster$IndexTable)
    terms.p$terms <- row.names(terms.p)
    terms.p.long <- melt(terms.p)
    g.N.long$p.val <- -log10(terms.p.long$value)
    #pdf(file=paste0(tag, "_enrich_", categoryType, "_lev3_heatmap.pdf"), height=9, width=10)
    #
    #nonZero.x <- subset(g.N.long, value>0)
    #g.out <- ggplot(nonZero.x, aes(variable, terms)) +
    #  geom_point(aes(colour=value, size=p.val)) +
    #  scale_color_gradient(low="black", high="steelblue") +
    #  xlab("Groups") + ylab(categoryType) +
    #  labs(colour="Gene numbers", size="-Log10(P-val)") +
    #  theme(panel.background=element_blank()) +
    #  theme(axis.text.x = element_text(colour = "black", angle=45, hjust=1),
    #    axis.text.y = element_text(colour = "black"))
    #print(g.out)
    #dev.off()

    write.table(g.N.long,
      file=paste0(tag, "_enrich_", categoryType, "_lev3_long_form.txt"),
      sep="\t", quote=FALSE, row.names=FALSE)
}

genes.list <- list()
genes.list <- lapply(config, function(x) {y <- readGeneList(x, e2s); y})

categoryTypes <- c("GO","KEGG","DOLITE","reactome.path")
# categoryTypes <- c("reactome.path")

# lapply(categoryTypes, gAnAnalyze, tag=tag, genes.list=genes.list,
#   species.anno=species.anno, pvalueT=pvalueT)

# categoryTypes <- c("GO", "GO.BP")

lapply(categoryTypes, gAnAnalyze.lev3, tag=tag, genes.list=genes.list,
  species.anno=species.anno, pvalueT=pvalueT)
