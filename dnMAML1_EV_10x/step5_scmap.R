library(Seurat)
library(SingleCellExperiment)
library(scmap)
load("All_Exp0.Rdata") # seurat object is called IMR90
# loading 10x data from senescence and GFP co-culture
load("imr90sen_growing_487cell_seurat.Robj")  # seurat object is called imr90sen_growing


met<-IMR90@meta.data
met$ID<-rep("mVenus:dnMAML1",571)
met$ID[met$ID_Group.celltype=="dnMAML1_RIS"]<-"RIS from dnMAML1 co-culture"
met$ID[met$ID_Group.celltype=="EV_mVenus"]<-"mVenus:EV"
met$ID[met$ID_Group.celltype=="EV_RIS"]<-"RIS from EV co-culture"
IMR90 <- AddMetaData(object = IMR90, metadata = met)   

# setting up dnMAML1 data
sce_dn <- SingleCellExperiment(assays = list(logcounts = as.matrix(IMR90@data)), colData = IMR90@meta.data$ID)
rowData(sce_dn)$feature_symbol <- rownames(sce_dn)

# remove features with duplicated names
sce_dn <- sce_dn[!duplicated(rownames(sce_dn)), ]
sce_dn
sce_dn <- selectFeatures(sce_dn, suppress_plot = FALSE)
table(rowData(sce_dn)$scmap_features)

#####
#  setting up co-culture data
sce_sen_growing <- SingleCellExperiment(assays = list(logcounts = as.matrix(imr90sen_growing@data)), colData = imr90sen_growing@meta.data$cluster)
rowData(sce_sen_growing)$feature_symbol <- rownames(sce_sen_growing)

# remove features with duplicated names
sce_sen_growing <- sce_sen_growing[!duplicated(rownames(sce_sen_growing)), ]
sce_sen_growing
sce_sen_growing <- selectFeatures(sce_sen_growing, suppress_plot = FALSE)
table(rowData(sce_sen_growing)$scmap_features)


####
#projection


sce_sen_growing <- indexCluster(sce_sen_growing,cluster_col="X")
# projection of dnmaml1 experiment to gfp+ris experiment
scmapCluster_results <- scmapCluster(
  projection = sce_dn, 
  index_list = list(
    imr90sen_growing = metadata(sce_sen_growing)$scmap_cluster_index
  )
)

head(scmapCluster_results$scmap_cluster_labs)
head(scmapCluster_results$scmap_cluster_siml)
head(scmapCluster_results$combined_labs)

plot(
  getSankey(
    colData(sce_dn)$X, 
    scmapCluster_results$scmap_cluster_labs[,'imr90sen_growing'],
    plot_height = 400
  )
)

###projecting the color to tsne plot
met<-IMR90@meta.data
met$project<-scmapCluster_results$scmap_cluster_labs
IMR90 <- AddMetaData(object = IMR90, metadata = met) 
p3<-TSNEPlot(object = IMR90,group.by="project",colors.use=c("grey","tan3","plum","gold"),pt.size=2)
p1<-TSNEPlot(object = IMR90,do.label = T,pt.size=2,group.by="ID_Group.celltype",colors.use = c("dodgerblue1","darkcyan","firebrick","aquamarine1"))
p2<-TSNEPlot(object = IMR90,do.label = T,pt.size=2,colors.use = c("burlywood1","chocolate1","chartreuse","darkslategray1"))
plot_grid(p1, p2,p3)

table(met$ID_Group.celltype,met$project)


slices <- c(247,80,4,3) 
lbls <- c("Growing","Secondary senescence","OIS","unassigned")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("grey","tan3","plum","gold"),
    main="mVenus:dnMAML1")

slices <- c(0,0,3,0) 
lbls <- c("Growing","Secondary senescence","OIS","unassigned")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("grey","tan3","plum","gold"),
    main="RIS from dnMAML1 co-culture")


slices <- c(128,81,2,8) 
lbls <- c("Growing","Secondary senescence","OIS","unassigned")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("grey","tan3","plum","gold"),
    main="mVenus:EV")

chisq.test(matrix(c(247,80,128,81),ncol=2))


slices <- c(0,0,15,0) 
lbls <- c("Growing","Secondary senescence","OIS","unassigned")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("grey","tan3","plum","gold"),
    main="RIS from EV co-culture")




# separate clusters of (1) secondary senescence into dnmaml1 and EV, (2) sen2, (3) growing
met<-IMR90@meta.data
met$final_cluster<-rep("Growing",571)
met$final_cluster[met$project=="OIS"]<-"OIS"
met$final_cluster[(met$project=="Secondary senescence") & (met$ID=="mVenus:EV")]<- "Secondary senescence, mVenus:EV"
met$final_cluster[(met$project=="Secondary senescence") & (met$ID=="mVenus:dnMAML1")]<- "Secondary senescence, mVenus:dnMAML1"
met$final_cluster[met$project=="unassigned"]<-"Unassigned"
IMR90 <- AddMetaData(object = IMR90, metadata = met) 

TSNEPlot(object = IMR90,group.by="final_cluster",colors.use=c("grey","plum","orange","red4","black"),pt.size=2)
## SCDE
Sen2<-data.frame(as.matrix(IMR90@raw.data))[,which(colnames(data.frame(as.matrix(IMR90@raw.data)))%in% rownames(met[met$final_cluster=="OIS",]))]
Growing<-data.frame(as.matrix(IMR90@raw.data))[,which(colnames(data.frame(as.matrix(IMR90@raw.data)))%in% rownames(met[met$final_cluster=="Growing",]))]
Sen1EV<-data.frame(as.matrix(IMR90@raw.data))[,which(colnames(data.frame(as.matrix(IMR90@raw.data)))%in% rownames(met[met$final_cluster=="Secondary senescence, mVenus:EV",]))]
Sen1dnMAML1<-data.frame(as.matrix(IMR90@raw.data))[,which(colnames(data.frame(as.matrix(IMR90@raw.data)))%in% rownames(met[met$final_cluster=="Secondary senescence, mVenus:dnMAML1",]))]


colnames(Sen2) <- paste("Sen2", colnames(Sen2), sep = "_")
colnames(Growing) <- paste("Growing", colnames(Growing), sep = "_")
colnames(Sen1EV) <- paste("Sen1EV", colnames(Sen1EV), sep = "_")
colnames(Sen1dnMAML1) <- paste("Sen1dnMAML1", colnames(Sen1dnMAML1), sep = "_")

Sen2_Growing<-cbind(Sen2,Growing)
Sen1EV_Growing<-cbind(Sen1EV,Growing)
Sen1dnMAML1_Growing<-cbind(Sen1dnMAML1,Growing)
Sen1dnMAML1_Sen1EV<-cbind(Sen1dnMAML1,Sen1EV)
Sen1EV_Sen2<-cbind(Sen1EV,Sen2)
Sen1dnMAML1_Sen2<-cbind(Sen1dnMAML1,Sen2)

write.table(Sen2_Growing,file=paste("./SCDE_all_exp0_projection/Sen2_Growing_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen1EV_Growing,file=paste("./SCDE_all_exp0_projection/Sen1EV_Growing_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen1dnMAML1_Growing,file=paste("./SCDE_all_exp0_projection/Sen1dnMAML1_Growing_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen1dnMAML1_Sen1EV,file=paste("./SCDE_all_exp0_projection/Sen1dnMAML1_Sen1EV_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen1EV_Sen2,file=paste("./SCDE_all_exp0_projection/Sen1EV_Sen2_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen1dnMAML1_Sen2,file=paste("./SCDE_all_exp0_projection/Sen1dnMAML1_Sen2_expr.txt",sep=""),sep="\t",quote=F)


##################################

library(scde)
# scde function
#expr matrix (cells on columns and genes on rows, cells on columns have to start with Cluster[0-9]_cellname)
runscde<-function(cluster1, cluster2){
  # SCDE
  library(scde)
  cluster1_2_expr<-read.table(paste("./SCDE_all_exp0_projection/",cluster1,"_",cluster2,"_expr.txt",sep=""),sep="\t",header=T,row.names = 1)
  
  con_cluster1_2 <- factor(gsub("_.*", "", colnames(cluster1_2_expr)), levels  =  c(paste(cluster1),paste(cluster2))) 
  # DE change is Cluster 1/2. Positive fold change = higher expression in Cluster 1
  names(con_cluster1_2) <- colnames(cluster1_2_expr)  
  table(con_cluster1_2)
  
  
  cluster1_2_expr<-apply(cluster1_2_expr,2,function(x) {storage.mode(x) <- 'integer'; x})

  o.ifm_cluster1_2split <- scde.error.models(counts = cluster1_2_expr, groups = con_cluster1_2, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1,min.count.threshold=1)
  
  
  
  # filter out cells that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.cells_cluster1_2 <- o.ifm_cluster1_2split$corr.a > 0
  table(valid.cells_cluster1_2)
  o.ifm_cluster1_2split <- o.ifm_cluster1_2split[valid.cells_cluster1_2, ]
  
  # estimate gene expression prior
  o.prior_cluster1_2 <- scde.expression.prior(models = o.ifm_cluster1_2split, counts = cluster1_2_expr, length.out = 400, show.plot = FALSE)
  
  
  
  
  #Testing for differential expression cluster1_2 split
  
  
  # define two groups_cluster1_2 of cells
  groups_cluster1_2 <- factor(gsub("_.*", "", rownames(o.ifm_cluster1_2split)), levels  =  c(paste(cluster1),paste(cluster2)))
  names(groups_cluster1_2) <- row.names(o.ifm_cluster1_2split)
  # run differential expression tests on all genes.
  ediff_cluster1_2 <- scde.expression.difference(o.ifm_cluster1_2split, cluster1_2_expr, o.prior_cluster1_2, groups  =  groups_cluster1_2, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
  
  # convert z score to p value

  p.values_cluster1_2 <- 2*pnorm(abs(ediff_cluster1_2$Z),lower.tail=F) # 2-tailed p-value
  p.values_cluster1_2.adj <- 2*pnorm(abs(ediff_cluster1_2$cZ),lower.tail=F) # Adjusted to control for FDR
  significant.genes_cluster1_2 <- which(p.values_cluster1_2.adj<1.1)
  length(significant.genes_cluster1_2)
  
  
  ord <- order(p.values_cluster1_2.adj[significant.genes_cluster1_2]) # order by p-value
  de_cluster1_2 <- cbind(ediff_cluster1_2[significant.genes_cluster1_2,1:3],p.values_cluster1_2.adj[significant.genes_cluster1_2])[ord,]
  colnames(de_cluster1_2) <- c("Lowerbound","log2FC","Upperbound","pvalue")
  write.table(de_cluster1_2,file=paste("./SCDE_all_exp0_projection/SCDE_",cluster1, "_",cluster2, ".txt",sep=""),quote=F,col.names=NA,sep="\t")
  
  de_cluster1_2<-data.frame(gene=rownames(de_cluster1_2),log2FC=de_cluster1_2$log2FC)
  write.table(de_cluster1_2,file=paste("./SCDE_all_exp0_projection/SCDE_",cluster1, "_",cluster2, ".rnk",sep=""),row.names=F,quote=F,col.names=NA,sep="\t")

}




runscde("Sen2","Growing")
runscde("Sen1EV","Growing")
runscde("Sen1dnMAML1","Growing")
runscde("Sen1dnMAML1","Sen1EV")
runscde("Sen1EV","Sen2")
runscde("Sen1dnMAML1","Sen2")


# the *rnk files from SCDE is used in GSEA, using hallmark.gmt as the gene sets

# heatmap for figure 3f
DoHeatmap(IMR90, cells.use = rownames(met[grep("Secondary",met$final_cluster),]), genes.use =c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A3","COL6A1","FGF7"),slim.col.label = TRUE,disp.min = -3,
          disp.max =3,col.low="blue",col.mid="white",col.high="red",group.spacing = 0.5,group.by = "final_cluster")

