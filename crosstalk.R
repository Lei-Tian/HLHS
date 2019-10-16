#For each gene, calculate the percetage of cells have expression>0 in Endocarduim and CM
myfunc=function(x){
	return(length(which(x>0))/length(x))
}
###############################################################
load("../Seurat/d83_LVp_LVn/all.Robj")

expr=all@data
ident=all@ident
###############################################################
CM=c(1,13)
Endocarduim=c(3,16)
clusters=list(CM,Endocarduim)
mynames=c("CM","Endocarduim")
for(i in c(1:2)){
	print(i)
	subCells=names(ident)[which(ident %in% clusters[[i]])]
	subExpr=expr[,subCells]
	pro=apply(subExpr,1,myfunc)
	x=names(pro)[which(pro>0.2)]
	write.table(x,paste("d83_LV/",mynames[i],".genes.0.2.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
}

ligandReceptor=read.table("pairs.txt",sep="\t")
res=data.frame(gene=c(names(cm),names(pro)),proportion=c(cm,pro),
               group=c(rep("CM",length(cm)),rep("EC",length(pro))))
res$gene_type="Other"
ligand=which(res$gene %in% ligandReceptor[,1])
res$gene_type[ligand]="Ligand"
receptor=which(res$gene %in% ligandReceptor[,2])
res$gene_type[receptor]="Receptor"
res$class=paste(res$group,res$gene_type,sep="_")
ggplot(res[which(res$gene_type!="Other"),],aes(x=proportion,colour=class))+geom_density()
