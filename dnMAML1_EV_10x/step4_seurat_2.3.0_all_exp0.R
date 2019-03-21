#Run step3_mVenus_puro_neo_exp_0.R first, some of the variables from that script are used in this script
library(Seurat)
library(dplyr)
library(Matrix)

# Load the dnMAML1 dataset
dnMAML1.data <- Read10X(data.dir ="dnMAML1_mVenus\\outs\\filtered_gene_bc_matrices\\hg19_mVenus_puro_neo\\")


x = as.matrix(x = dnMAML1.data)

#dnMAML1_RIS_mVenus<-dnMAML1_RIS_mVenus[grep("mVenus",dnMAML1_RIS_mVenus$Con1),]
x.subset<-x[,colnames(x)%in% rownames(dnMAML1_RIS_mVenus)]


x.subset<-x.subset[,order(colnames(x.subset))]

dnMAML1_RIS_mVenus<-dnMAML1_RIS_mVenus[rownames(dnMAML1_RIS_mVenus) %in% colnames(x.subset),]

dnMAML1_RIS_mVenus<-dnMAML1_RIS_mVenus[order(rownames(dnMAML1_RIS_mVenus)),]
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
dnMAML1 <- CreateSeuratObject(raw.data = x.subset,min.cells=3,min.genes=200,project="10X_dnMAML1")




# Load the EV dataset
EV.data <- Read10X(data.dir ="EV_mVenus\\outs\\filtered_gene_bc_matrices\\hg19_mVenus_puro_neo\\")

x = as.matrix(x = EV.data)
x.subset<-x[,colnames(x)%in% rownames(EV_RIS_mVenus)]



x.subset<-x.subset[,order(colnames(x.subset))]

EV_RIS_mVenus<-EV_RIS_mVenus[rownames(EV_RIS_mVenus) %in% colnames(x.subset),]

EV_RIS_mVenus<-EV_RIS_mVenus[order(rownames(EV_RIS_mVenus)),]


EV <- CreateSeuratObject(raw.data = x.subset,min.cells=3,min.genes=200,project="10X_EV")





IMR90 <- MergeSeurat(object1 = EV, object2 = dnMAML1, add.cell.id1 = "EV", add.cell.id2 = "dnMAML1", project = "IMR90")


mito.genes <- grep(pattern = "^MT-", x = rownames(x = IMR90@data), value = TRUE)

percent.mito <- Matrix::colSums(IMR90@raw.data[mito.genes, ])/Matrix::colSums(IMR90@raw.data)

IMR90 <- AddMetaData(object = IMR90, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = IMR90, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = IMR90, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = IMR90, gene1 = "nUMI", gene2 = "nGene")



IMR90 <- FilterCells(object = IMR90, subset.names = c("nGene", "percent.mito"), 
                     low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.1))


# normalize
IMR90 <- NormalizeData(object = IMR90, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# find variable gene
IMR90 <- FindVariableGenes(object = IMR90, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

########################
IMR90 <- ScaleData(object = IMR90, vars.to.regress = c("nUMI", "percent.mito"))


IMR90 <- RunPCA(object = IMR90, pc.genes = IMR90@var.genes, do.print = TRUE, pcs.print = 1:5, 
                genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = IMR90, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = IMR90, pcs.use = 1:2)

IMR90 <- AddMetaData(object = IMR90, metadata = passed)



PCAPlot(object = IMR90, dim.1 = 1, dim.2 = 2,group.by="ID_Group.celltype")

IMR90 <- ProjectPCA(object = IMR90, do.print = FALSE)
PCHeatmap(object = IMR90, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = IMR90, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


#Determine statistically significant principal components

IMR90 <- JackStraw(object = IMR90, num.replicate = 100)



###########################################
JackStrawPlot(IMR90, PCs = 1:19)
PCElbowPlot(IMR90, num.pc = 40)
# The clustering from Seurat is not used. 



IMR90 <- FindClusters(object = IMR90, reduction.type = "pca", dims.use = 1:7, 
                      resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc = T)

# TSNE
IMR90 <- RunTSNE(object = IMR90, dims.use = 1:7, do.fast = TRUE)




save(IMR90,file="All_Exp0.Rdata")


