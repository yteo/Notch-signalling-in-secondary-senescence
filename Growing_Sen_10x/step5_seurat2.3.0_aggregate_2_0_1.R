library(Seurat)
library(dplyr)
library(Matrix)
# seurat 2.0.1
# Load the imr90sen_growing dataset
imr90sen_growing.data <- Read10X(data.dir ="IMR90_sen_growing_libraries_2_0_1\\outs\\filtered_gene_bc_matrices_mex\\hg19_neo_gfp_puro\\")

x = as.matrix(x = imr90sen_growing.data)
passed<-read.table("passedcell_GFP_RIS_stringentGFP.txt",header=T)
#passed<-passed[grep("sen1|growing",passed$Con2),]
names.use <- colnames(x)[(colnames(x) %in% rownames(passed))]
x.subset <- x[, names.use]
x.subset<-x.subset[,order(colnames(x.subset))]
passed<-passed[order(rownames(passed)),]

imr90sen_growing <- CreateSeuratObject(raw.data = x.subset,min.cells=3,min.genes=200,project="10X_IMR90sengrowing")


#####################
# basic QC and selecting cells for further analysis

mito.genes <- grep(pattern = "^MT-", x = rownames(x = imr90sen_growing@data), value = TRUE)
percent.mito <- colSums(imr90sen_growing@raw.data[mito.genes, ])/colSums(imr90sen_growing@raw.data)


imr90sen_growing <- AddMetaData(object = imr90sen_growing, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = imr90sen_growing, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = imr90sen_growing, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = imr90sen_growing, gene1 = "nUMI", gene2 = "nGene")


imr90sen_growing <- FilterCells(object = imr90sen_growing, subset.names = c("nGene", "percent.mito"), 
                                low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.15))
#####################


imr90sen_growing <- AddMetaData(object = imr90sen_growing, metadata = passed, col.name = c("Con1","Con2"))

# normalize
imr90sen_growing <- NormalizeData(object = imr90sen_growing, normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
# find variable gene
imr90sen_growing <- FindVariableGenes(object = imr90sen_growing, mean.function = ExpMean, dispersion.function = LogVMR, 
                                      x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)


cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

imr90sen_growing <- CellCycleScoring(object = imr90sen_growing, s.genes = s.genes, g2m.genes = g2m.genes, 
                                     set.ident = TRUE)


imr90sen_growing@meta.data$CC.Difference <- imr90sen_growing@meta.data$S.Score - imr90sen_growing@meta.data$G2M.Score

# scaling and regress out effects
imr90sen_growing <- ScaleData(object = imr90sen_growing, vars.to.regress = c("nUMI", "percent.mito","CC.Difference"))
########################

#Perform linear dimensional reduction
imr90sen_growing <- RunPCA(object = imr90sen_growing, pc.genes = imr90sen_growing@var.genes, do.print = TRUE, pcs.print = 1:5, 
                           genes.print = 5)
# Examine and visualize PCA results a few different ways
PrintPCA(object = imr90sen_growing, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = imr90sen_growing, pcs.use = 1:2)
PCAPlot(object = imr90sen_growing, dim.1 = 1, dim.2 = 2,group.by="Con2")


imr90sen_growing <- ProjectPCA(object = imr90sen_growing, do.print = FALSE)
PCHeatmap(object = imr90sen_growing, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = imr90sen_growing, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)



imr90sen_growing <- JackStraw(object = imr90sen_growing, num.replicate = 100, do.print = FALSE)




###########################################
JackStrawPlot(imr90sen_growing, PCs = 1:19)
PCElbowPlot(imr90sen_growing, num.pc = 40)
# cluster cellls


# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
imr90sen_growing <- FindClusters(object = imr90sen_growing, reduction.type = "pca", dims.use = 1:15, 
                                 resolution = 0.4, print.output = 0, save.SNN = TRUE,force.recalc = T)

# TSNE
imr90sen_growing <- RunTSNE(object = imr90sen_growing, dims.use = 1:15, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = imr90sen_growing,group.by="Con3",colors.use=c("black","skyblue","brown1"),pt.size=2)


################
library(sparcl)
#extract tsne 1 and 2 and do own cluster
set.seed(656)
tsne_imr90<-data.frame(imr90sen_growing@dr$tsne@cell.embeddings)


tsne_imr90_original=tsne_imr90
## Creating k-means clustering model and assigning to the data used to create the tsne
km.perm=KMeansSparseCluster.permute(tsne_imr90,K=3,nperms=100)
fit_cluster_kmeans<- KMeansSparseCluster(tsne_imr90,K=3,wbounds=km.perm$bestw)

tsne_imr90_original$kmeanscluster = factor(fit_cluster_kmeans[[1]]$Cs)


ggplot(tsne_imr90_original, aes_string(x="tSNE_1", y="tSNE_2", color="kmeanscluster")) +
  geom_point(size=2) +scale_color_manual(values=c("burlywood3","grey","purple","green","cyan","darkgoldenrod"))

imr90sen_growing <- AddMetaData(object = imr90sen_growing, metadata = tsne_imr90_original)    
met<-imr90sen_growing@meta.data

table(met$Con3,met$kmeanscluster)

TSNEPlot(object = imr90sen_growing,group.by="kmeanscluster",colors.use=c("plum","grey","tan3"),pt.size=2)
save(imr90sen_growing, file = "imr90sen_growing_487cell_seurat.Robj")



chisq.test(table(met$Con3,met$kmeanscluster)[,c(1,3)][c(1,3),]) 




slices <- c(28,3) 
lbls <- c("Secondary senescence", "OIS")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("yellow","plum"),
    main="Day7 GFP")


slices <- c(50,177) 
lbls <- c("Secondary senescence", "OIS")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=c("yellow","plum"),
    main="Day7 OIS")
######################


# use the boxplot function at the end of this script
Boxplot(imr90sen_growing,c("CDKN1A","CDKN2B","IL1B","IL6","IL8"),cols.use = c("white","white","white"))
Boxplot(imr90sen_growing,c("COL1A1","COL3A1","COL5A2","TGFB1I1","CEBPB","CTGF"),cols.use = c("white","white","white"),group.by="kmeanscluster")

######################



met<-imr90sen_growing@meta.data



#####################
# cluster 1 consists of majority of  OIS cells
# cluster 2 consists of growing cells
# cluster 3 consists of secondary senescence cells

met<-imr90sen_growing@meta.data
met_R<-met[((met$kmeanscluster=="1")),]
met_N<-met[((met$kmeanscluster=="3")),]
met_G<-met[((met$kmeanscluster=="2")),]
met_R_N<-met[((met$kmeanscluster=="3")|(met$kmeanscluster=="1")),]


Cluster1R<-data.frame(as.matrix(imr90sen_growing@raw.data))[,which(colnames(data.frame(as.matrix(imr90sen_growing@raw.data)))%in% gsub("-",".",rownames(met_R)))]
Growing<-data.frame(as.matrix(imr90sen_growing@raw.data))[,which(colnames(data.frame(as.matrix(imr90sen_growing@raw.data)))%in% gsub("-",".",rownames(met_G)))]
Cluster3N<-data.frame(as.matrix(imr90sen_growing@raw.data))[,which(colnames(data.frame(as.matrix(imr90sen_growing@raw.data)))%in% gsub("-",".",rownames(met_N)))]

Cluster1R_3N<-data.frame(as.matrix(imr90sen_growing@raw.data))[,which(colnames(data.frame(as.matrix(imr90sen_growing@raw.data)))%in% gsub("-",".",rownames(met_R_N)))]

colnames(Cluster1R) <- paste("Cluster1R", colnames(Cluster1R), sep = "_")
colnames(Growing) <- paste("Growing", colnames(Growing), sep = "_")
colnames(Cluster3N) <- paste("Cluster3N", colnames(Cluster3N), sep = "_")
colnames(Cluster1R_3N) <- paste("Sen", colnames(Cluster1R_3N), sep = "_")

Cluster1R_Growing<-cbind(Cluster1R,Growing)
Cluster3N_Growing<-cbind(Cluster3N,Growing)
Cluster3N_Cluster1R<-cbind(Cluster3N,Cluster1R)
Sen_Growing<-cbind(Cluster1R_3N,Growing)



write.table(Cluster1R_Growing,file=paste("Cluster1R_Growing_expr.txt",sep=""),sep="\t",quote=F)
write.table(Cluster3N_Growing,file=paste("Cluster3N_Growing_expr.txt",sep=""),sep="\t",quote=F)
write.table(Cluster3N_Cluster1R,file=paste("Cluster3N_Cluster1R_expr.txt",sep=""),sep="\t",quote=F)
write.table(Sen_Growing,file=paste("Sen_Growing_expr.txt",sep=""),sep="\t",quote=F)

####################################


##################################

library(scde)
# scde function
#expr matrix (cells on columns and genes on rows, cells on columns have to start with Cluster[0-9]_cellname)
runscde<-function(cluster1, cluster2){
  # SCDE
  library(scde)
  cluster1_2_expr<-read.table(paste("",cluster1,"_",cluster2,"_expr.txt",sep=""),sep="\t",header=T,row.names = 1)
  
  con_cluster1_2 <- factor(gsub("_.*", "", colnames(cluster1_2_expr)), levels  =  c(paste(cluster1),paste(cluster2))) 
  # DE change is Cluster 1/2. Positive fold change = higher expression in Cluster 1
  names(con_cluster1_2) <- colnames(cluster1_2_expr)  
  table(con_cluster1_2)
  
  
  cluster1_2_expr<-apply(cluster1_2_expr,2,function(x) {storage.mode(x) <- 'integer'; x})
  #Can I use SCDE/PAGODA with UMI data?
  #A: Yes. In the error building step, you will need to reduce min.count.threshold to 1.
  #write.table(cluster1_2_expr,file="cluster1_2_expr.txt",sep="\t",quote=F)
  o.ifm_cluster1_2split <- scde.error.models(counts = cluster1_2_expr, groups = con_cluster1_2, n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1,min.count.threshold=1)
  
  
  
  # filter out cells that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.cells_cluster1_2 <- o.ifm_cluster1_2split$corr.a > 0
  table(valid.cells_cluster1_2)
  o.ifm_cluster1_2split <- o.ifm_cluster1_2split[valid.cells_cluster1_2, ]
  
  # estimate gene expression prior
  o.prior_cluster1_2 <- scde.expression.prior(models = o.ifm_cluster1_2split, counts = cluster1_2_expr, length.out = 400, show.plot = FALSE)
  
  
  
  
  #Testing for differential expression cluster1_2 split
  
  
  # define two groups_cluster1_2 of cells
  groups_cluster1_2 <- factor(gsub("_.*", "", rownames(o.ifm_cluster1_2split)), levels  =  c(paste(cluster1),paste(Cluster2)))
  names(groups_cluster1_2) <- row.names(o.ifm_cluster1_2split)
  # run differential expression tests on all genes.
  ediff_cluster1_2 <- scde.expression.difference(o.ifm_cluster1_2split, cluster1_2_expr, o.prior_cluster1_2, groups  =  groups_cluster1_2, n.randomizations  =  100, n.cores  =  4, verbose  =  1)
  
  # convert z score to p value
  
  #http://hms-dbmi.github.io/scw/differential-expression.html
  p.values_cluster1_2 <- 2*pnorm(abs(ediff_cluster1_2$Z),lower.tail=F) # 2-tailed p-value
  p.values_cluster1_2.adj <- 2*pnorm(abs(ediff_cluster1_2$cZ),lower.tail=F) # Adjusted to control for FDR
  significant.genes_cluster1_2 <- which(p.values_cluster1_2.adj<1.1)
  length(significant.genes_cluster1_2)
  
  
  ord <- order(p.values_cluster1_2.adj[significant.genes_cluster1_2]) # order by p-value
  de_cluster1_2 <- cbind(ediff_cluster1_2[significant.genes_cluster1_2,1:3],p.values_cluster1_2.adj[significant.genes_cluster1_2])[ord,]
  colnames(de_cluster1_2) <- c("Lowerbound","log2FC","Upperbound","pvalue")
  write.table(de_cluster1_2,file=paste("SCDE_",cluster1, "_",Cluster2, ".txt",sep=""),quote=F,col.names=NA,sep="\t")
  
  de_cluster1_2<-data.frame(gene=rownames(de_cluster1_2),log2FC=de_cluster1_2$log2FC)
  write.table(de_cluster1_2,file=paste("SCDE_",cluster1, "_",Cluster2, ".rnk",sep=""),row.names=F,quote=F,sep="\t")
  
  #DE_Cluster0_Cluster3<-data.frame(gene=rownames(de_cluster1_2),log2FC=de_cluster1_2$log2FC)
  
}




runscde("Cluster1R","Growing")
runscde("Cluster3N","Growing")
runscde("Cluster3N","Cluster1R")
runscde("Sen","Growing")

############## modified functions to do boxplot instead of vlnplot in seurat

SingleVlnPlot2<-function (feature, data, cell.ident, do.sort, y.max, size.x.use, 
                          size.y.use, size.title.use, adjust.use, point.size.use, cols.use, 
                          gene.names, y.log, x.lab.rot, y.lab.rot, legend.position, 
                          remove.legend) 
{
  set.seed(seed = 42)
  data$ident <- cell.ident
  if (do.sort) {
    data$ident <- factor(x = data$ident, levels = names(x = rev(x = sort(x = tapply(X = data[, 
                                                                                             feature], INDEX = data$ident, FUN = mean)))))
  }
  if (y.log) {
    noise <- rnorm(n = length(x = data[, feature]))/200
    data[, feature] <- data[, feature] + 1
  }
  else {
    noise <- rnorm(n = length(x = data[, feature]))/1e+05
  }
  data[, feature] <- data[, feature] + noise
  y.max <- SetIfNull(x = y.max, default = max(data[, feature]))
  data$ident <- factor(x = data$ident, levels =c(2,1,3))
  plot <- ggplot(data = data, mapping = aes(x = factor(x = ident), 
                                            y = eval(expr = parse(text = feature)))) + geom_boxplot(outlier.shape=NA ,adjust = adjust.use, trim = TRUE, mapping = aes(fill = factor(x = ident))) + 
    theme(legend.position = legend.position, axis.title.x = element_text(face = "bold", 
                                                                         colour = "#990000", size = size.x.use), axis.title.y = element_text(face = "bold", 
                                                                                                                                             colour = "#990000", size = size.y.use)) + guides(fill = guide_legend(title = NULL)) + 
    # geom_jitter(height = 0, size = point.size.use) + xlab("Identity") + 
    NoGrid() + ggtitle(feature) + theme(plot.title = element_text(size = size.title.use, 
                                                                  face = "bold"))
  if (y.log) {
    plot <- plot + scale_y_log10()
  }
  else {
    plot <- plot + ylim(min(data[, feature]), y.max)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + ylab(label = "Log Expression level")
    }
    else {
      plot <- plot + ylab(label = "Expression level")
    }
  }
  else {
    plot <- plot + ylab(label = "")
  }
  if (!is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, 
                                                    vjust = 0.5))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}
NoGrid<-function (...) 
{
  no.grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   ...)
  return(no.grid)
}
Boxplot<-function (object, features.plot, ident.include = NULL, nCol = NULL, 
                   do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16, 
                   size.y.use = 16, size.title.use = 20, adjust.use = 1, point.size.use = 1, 
                   cols.use = NULL, group.by = NULL, y.log = FALSE, x.lab.rot = FALSE, 
                   y.lab.rot = FALSE, legend.position = "right", single.legend = TRUE, 
                   remove.legend = FALSE, do.return = FALSE, return.plotlist = FALSE, 
                   ...) 
{
  if (is.null(x = nCol)) {
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
    else {
      nCol <- min(length(x = features.plot), 3)
    }
  }
  data.use <- data.frame(FetchData(object = object, vars.all = features.plot, 
                                   ...))
  if (is.null(x = ident.include)) {
    cells.to.include <- object@cell.names
  }
  else {
    cells.to.include <- WhichCells(object = object, ident = ident.include)
  }
  data.use <- data.use[cells.to.include, , drop = FALSE]
  if (!is.null(x = group.by)) {
    ident.use <- as.factor(x = FetchData(object = object, 
                                         vars.all = group.by)[cells.to.include, 1])
  }
  else {
    ident.use <- object@ident[cells.to.include]
  }
  gene.names <- colnames(x = data.use)[colnames(x = data.use) %in% 
                                         rownames(x = object@data)]
  if (single.legend) {
    remove.legend <- TRUE
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(X = features.plot, FUN = function(x) {
    return(SingleVlnPlot2(feature = x, data = data.use[, x, 
                                                       drop = FALSE], cell.ident = ident.use, do.sort = do.sort, 
                          y.max = y.max, size.x.use = size.x.use, size.y.use = size.y.use, 
                          size.title.use = size.title.use, adjust.use = adjust.use, 
                          point.size.use = point.size.use, cols.use = cols.use, 
                          gene.names = gene.names, y.log = y.log, x.lab.rot = x.lab.rot, 
                          y.lab.rot = y.lab.rot, legend.position = legend.position, 
                          remove.legend = remove.legend))
  })
  if (length(x = features.plot) > 1) {
    plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
    if (single.legend && !remove.legend) {
      legend <- get_legend(plot = plots[[1]] + theme(legend.position = legend.position))
      if (legend.position == "bottom") {
        plots.combined <- plot_grid(plots.combined, legend, 
                                    ncol = 1, rel_heights = c(1, 0.2))
      }
      else if (legend.position == "right") {
        plots.combined <- plot_grid(plots.combined, legend, 
                                    rel_widths = c(3, 0.3))
      }
      else {
        warning("Shared legends must be at the bottom or right of the plot")
      }
    }
  }
  else {
    plots.combined <- plots[[1]]
  }
  if (do.return) {
    if (return.plotlist) {
      return(plots)
    }
    else {
      return(plots.combined)
    }
  }
  else {
    if (length(x = plots.combined) > 1) {
      plots.combined
    }
    else {
      invisible(x = lapply(X = plots.combined, FUN = print))
    }
  }
}
SetIfNull<-function (x, default) 
{
  if (is.null(x = x)) {
    return(default)
  }
  else {
    return(x)
  }
}
#################


# venn diagram 
a=paracrine
b=coculture secondary / coculture OIS
c=timecourse secondary / timecourse OIS

library(venneuler)
v <- venneuler(c(a=715, b=245,c=17 ,"a&b"=43,"a&c"=3,"b&c"=21,"a&b&c"=14))
plot(v)


