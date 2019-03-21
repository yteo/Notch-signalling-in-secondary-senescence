library(Seurat)
# reading in smart seq data and the branching of senescent cells
# senbotoom in the smart seq data is the secondary senescence branch
# sentop is the OIS branch
dat<-read.table("filtered_raw_IMR90.txt")
rownames(dat)<-make.names(dat$name,unique=T)
dat<-dat[,-1]
splitcells<-read.csv("split.csv",header=T)

bottom<-dat[,colnames(dat) %in% splitcells$Sen_17cells_individualbranch]
top<-dat[,colnames(dat) %in% splitcells$Sen_22cells_day24branch]
colnames(bottom)<-paste("SenBottom",colnames(bottom),sep="_")
colnames(top)<-paste("SenTop",colnames(top),sep="_")
smartseq.data<-cbind(bottom,top)
########
# 10x (reading 10x data from co-culture)
tenx.data <- Read10X(data.dir ="IMR90_sen_growing_libraries_2_0_1\\outs\\filtered_gene_bc_matrices_mex\\hg19_neo_gfp_puro\\")

x = as.matrix(x = tenx.data)
passed<-read.table("passedcell_GFP_RIS_stringentGFP.txt",header=T)

load("imr90sen_growing_487cell_seurat.Robj")

met<-imr90sen_growing@meta.data

met3<-met[rownames(met) %in% rownames(passed),]


met3$SenClust[met3$kmeanscluster==1]<-"OIS"
met3$SenClust[met3$kmeanscluster==3]<-"Secondary"

met3<-met3[,-1]
met3<-met3[,-1]
met3<-met3[,-1]
met3<-met3[,-1]
met3<-met3[which(met3$Con3!="Growing"),]
met3<-met3[complete.cases(met3), ]
names.use <- colnames(x)[(colnames(x) %in% rownames(met3))]
x.subset <- x[, names.use]
x.subset<-x.subset[,order(colnames(x.subset))]
met3<-met3[order(rownames(met3)),]




############################################


smartseq <- CreateSeuratObject(raw.data = smartseq.data)
smartseq <- NormalizeData(object = smartseq)
met<-smartseq@meta.data
mito.genes <- grep(pattern = "^MT-", x = rownames(x = smartseq.data), value = TRUE)

smartseq <- ScaleData(object = smartseq, vars.to.regress = c("nUMI"))
smartseq <- FindVariableGenes(object = smartseq, do.plot = FALSE)


##############################################
tenx <- CreateSeuratObject(raw.data = x.subset,min.cells=3,min.genes=200)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = tenx@data), value = TRUE)
percent.mito <- colSums(tenx@raw.data[mito.genes, ])/colSums(tenx@raw.data)


tenx <- AddMetaData(object = tenx, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = tenx, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = tenx, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = tenx, gene1 = "nUMI", gene2 = "nGene")



tenx <- FilterCells(object = tenx, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.15))
#####################


tenx <- AddMetaData(object = tenx, metadata = met3)

# normalize
tenx <- NormalizeData(object = tenx, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
# find variable gene
tenx <- FindVariableGenes(object = tenx, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
tenx <- ScaleData(object = tenx, vars.to.regress = c("nUMI", "percent.mito"))


##########################
# we will take the union of the top 50 variable genes in each dataset for alignment
hvg.smartseq <- rownames(x = head(x = smartseq@hvg.info, n = 50))
hvg.tenx <- rownames(x = head(x = tenx@hvg.info, n = 50))
hvg.union <- union(x = hvg.smartseq, y = hvg.tenx)

tenx@meta.data[, "protocol"] <- "10X"
smartseq@meta.data[, "protocol"] <- "smartseq"


sen <- RunCCA(object = tenx, object2 = smartseq, genes.use = hvg.union)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = sen , reduction.use = "cca", group.by = "protocol", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = sen , features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2)


PrintDim(object = sen, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
DimHeatmap(object = sen, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

DimHeatmap(object = sen, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

sen <- CalcVarExpRatio(object = sen, reduction.type = "pca", grouping.var = "protocol", 
                       dims.use = 1:2)

# We discard cells where the variance explained by CCA is <2-fold (ratio <0.5) compared to PCA
sen.all.save <- sen
sen <- SubsetData(object = sen, subset.name = "var.ratio.pca", accept.low = 0.5)

sen.discard <- SubsetData(object = sen.all.save, subset.name = "var.ratio.pca", 
                          accept.high = 0.5)
median(x = sen@meta.data[, "nGene"])

median(x = sen.discard@meta.data[, "nGene"])

sen <- AlignSubspace(object = sen, reduction.type = "cca", grouping.var = "protocol", 
                     dims.align = 1:2)
#Visualize the aligned CCA and perform integrated analysis

p1 <- VlnPlot(object = sen, features.plot = "ACC1", group.by = "protocol", 
              do.return = TRUE)
p2 <- VlnPlot(object = sen, features.plot = "ACC2", group.by = "protocol", 
              do.return = TRUE)
plot_grid(p1, p2)

sen <- RunTSNE(object = sen, reduction.use = "cca.aligned", dims.use = 1:2, 
               do.fast = TRUE)
sen <- FindClusters(object = sen, reduction.type = "cca.aligned", dims.use = 1:2,  resolution=0.3,
                    save.SNN = TRUE,force.recalc=T)
met<-sen@meta.data
met$Con4<-as.character(met$SenClust)
met$Con4[met$orig.ident=="SenBottom"]<-"SenBottom"
met$Con4[met$orig.ident=="SenTop"]<-"SenTop"
sen <- AddMetaData(object = sen, metadata = met)


TSNEPlot(object = sen, group.by = "Con4", do.return = TRUE, pt.size = 2,colors.use = c("blue","tan3","plum4","deeppink1"))


library(sparcl)
set.seed(656)
tsne_imr90<-data.frame(sen@dr$tsne@cell.embeddings)

# kmeans

tsne_imr90_original=tsne_imr90
## Creating k-means clustering model
km.perm=KMeansSparseCluster.permute(tsne_imr90,K=2,nperms=100)
fit_cluster_kmeans<- KMeansSparseCluster(tsne_imr90,K=2,wbounds=km.perm$bestw)

tsne_imr90_original$kmeanscluster = factor(fit_cluster_kmeans[[1]]$Cs)


ggplot(tsne_imr90_original, aes_string(x="tSNE_1", y="tSNE_2", color="kmeanscluster")) +
  geom_point(size=2) +scale_color_manual(values=c("burlywood3","grey","purple","green","cyan","darkgoldenrod"))


geom_point(size=2) +scale_color_manual(values=c("burlywood3","grey","purple","green","cyan","darkgoldenrod"))
sen <- AddMetaData(object = sen, metadata = tsne_imr90_original)    
met<-sen@meta.data

table(met$Con4,met$kmeanscluster)

TSNEPlot(object = sen,group.by="kmeanscluster",colors.use=c("plum","grey","tan3"),pt.size=2)



# secondary and OIS in this case are from co-culture experiment
# senescent top and bottom are from time-course experiment


slices <- c(3, 7,119,18) 
lbls <- c("Senescent bottom", "Secondary", "OIS","Senescent top")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("blue","cyan","darkgoldenrod3","deeppink1"),
    main="1")


slices <- c(12,70,33,2) 
lbls <- c("Senescent bottom", "Secondary", "OIS","Senescent top")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("blue","cyan","darkgoldenrod3","deeppink1"),
    main="2")



# chi-square
 #               1   2
 # SenBottom     3  12
 # SenTop       18   2


 #               1   2
 # Secondary   7  70
 # OIS 119  33

library(fifer)
chisq.post.hoc(table(met$Con4,met$kmeanscluster), test='chisq.test')
chisq.test(table(met$Con4,met$kmeanscluster)[c(1,4),])

chisq.test(table(met$Con4,met$kmeanscluster)[c(2,3),])
