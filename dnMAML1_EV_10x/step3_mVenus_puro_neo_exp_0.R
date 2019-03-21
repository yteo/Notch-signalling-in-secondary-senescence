library(Seurat)
library(dplyr)
library(Matrix)
# Load the dnMAML1 dataset
dnMAML1.data <- Read10X(data.dir ="dnMAML1_mVenus\\outs\\filtered_gene_bc_matrices\\hg19_mVenus_puro_neo\\")


dnMAML1 <- CreateSeuratObject(raw.data = dnMAML1.data,min.cells=3,min.genes=200,project="10X_dnMAML1")


#####################
# basic QC and selecting cells for further analysis


mito.genes <- grep(pattern = "^MT-", x = rownames(x = dnMAML1@data), value = TRUE)
percent.mito <- colSums(dnMAML1@raw.data[mito.genes, ])/colSums(dnMAML1@raw.data)


dnMAML1 <- AddMetaData(object = dnMAML1, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = dnMAML1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = dnMAML1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = dnMAML1, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts less than 2500 and percent mitochondria more than 0.1
dnMAML1 <- FilterCells(object = dnMAML1, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.1))


# normalize
dnMAML1 <- NormalizeData(object = dnMAML1, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

# find variable gene
dnMAML1 <- FindVariableGenes(object = dnMAML1, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
dnMAML1 <- ScaleData(object = dnMAML1, vars.to.regress = c("nUMI", "percent.mito"))

########################################
# finding dnMAML1_mVenus positive cells based on the expression of puro, neomycin or mVenus expression > 0
# ras+ based on freebayes (step2_ras_snp.sh)
dnMAML1_Ras<-data.frame(rownames=c("ACGAGGATCATGCAAC","AGATTGCGTACACCGC","ATCATGGTCGAATGCT","CAAGTTGGTGAGCGAT","CGAGCCATCCAAATGC","CGATGTATCTGGCGAC","CGCTGGACACGAGGTA","GACTGCGTCAATCACG","GATCGTACAGTCCTTC","GCACATAAGCTGAAAT"
),Con1=rep("dnMAML1_RIS",10),Con2=rep("RIS",10))  # identified using step2_ras_snp.sh
rownames(dnMAML1_Ras)<-dnMAML1_Ras$rownames
dnMAML1_Ras<-dnMAML1_Ras[,-1]
dnMAML1_mVenus<-data.frame(t(data.frame(as.matrix(dnMAML1@data[grep("^mVenus$|^pLPC_puro$",rownames(dnMAML1@data)),]))))
dnMAML1_mVenus_pos<-data.frame(dnMAML1_mVenus[which(dnMAML1_mVenus$mVenus>0|dnMAML1_mVenus$pLPC_puro > 0),])
# finding dnMAML1_neo positive cells
dnMAML1_neo<-data.frame(dnMAML1@data[grep("^Neomycin$",rownames(dnMAML1@data)),])
dnMAML1_neo<-data.frame(dnMAML1_neo[which(dnMAML1_neo$dnMAML1.data.grep...Neomycin....rownames.dnMAML1.data.....>0),,drop=F])
dnMAML1_neo_pos<-dnMAML1_neo
dnMAML1_Ras_pos<-dnMAML1_Ras
Reduce(intersect,list(rownames(dnMAML1_mVenus_pos),rownames(dnMAML1_neo),rownames(dnMAML1_Ras)))
intersect(rownames(dnMAML1_mVenus_pos),rownames(dnMAML1_neo))
intersect(rownames(dnMAML1_mVenus_pos),rownames(dnMAML1_Ras))
intersect(rownames(dnMAML1_neo),rownames(dnMAML1_Ras))

# mVenus only cell
dnMAML1_mVenus<-subset(dnMAML1_mVenus_pos, !(rownames(dnMAML1_mVenus_pos) %in% rownames(dnMAML1_neo)))
dnMAML1_mVenus<-subset(dnMAML1_mVenus, !(rownames(dnMAML1_mVenus) %in% rownames(dnMAML1_Ras)))

# of rasv12 only cell
dnMAML1_Ras<-subset(dnMAML1_Ras, !(rownames(dnMAML1_Ras) %in% rownames(dnMAML1_mVenus_pos)))
dnMAML1_Ras<-subset(dnMAML1_Ras, !(rownames(dnMAML1_Ras) %in% rownames(dnMAML1_neo_pos)))

dnMAML1_neo<-subset(dnMAML1_neo, !(rownames(dnMAML1_neo) %in% rownames(dnMAML1_mVenus_pos)))
# of neo only cells
dnMAML1_neo<-subset(dnMAML1_neo, !(rownames(dnMAML1_neo) %in% rownames(dnMAML1_Ras_pos)))

# neo and rasv12 only cells
dnMAML1_neo_ras<-subset(dnMAML1_neo_pos,rownames(dnMAML1_neo_pos) %in% rownames(dnMAML1_Ras_pos))
dnMAML1_neo_ras<-subset(dnMAML1_neo_ras,!(rownames(dnMAML1_neo_ras) %in% rownames(dnMAML1_mVenus_pos)))


dnMAML1_mVenus$Con1<-rep("dnMAML1_mVenus",dim(dnMAML1_mVenus)[1])
dnMAML1_neo$Con1<-rep("dnMAML1_RIS",dim(dnMAML1_neo)[1])
dnMAML1_mVenus$Con2<-rep("mVenus",dim(dnMAML1_mVenus)[1])
dnMAML1_neo$Con2<-rep("RIS",dim(dnMAML1_neo)[1])
dnMAML1_neo_ras$Con1<-rep("dnMAML1_RIS",dim(dnMAML1_neo_ras)[1])
dnMAML1_neo_ras$Con2<-rep("RIS",dim(dnMAML1_neo_ras)[1])

dnMAML1_neo_ras<-dnMAML1_neo_ras[,-1]
dnMAML1_mVenus<-dnMAML1_mVenus[,-c(1,2)]
dnMAML1_neo<-dnMAML1_neo[,-1]
dnMAML1_RIS_mVenus<-rbind(dnMAML1_mVenus,dnMAML1_neo,dnMAML1_Ras,dnMAML1_neo_ras)
###############################################################################################################
# EV
# Load the EV dataset
EV.data <- Read10X(data.dir ="EV_mVenus\\outs\\filtered_gene_bc_matrices\\hg19_mVenus_puro_neo\\")


EV <- CreateSeuratObject(raw.data = EV.data,min.cells=3,min.genes=200,project="10X_EV")


#####################
# basic QC and selecting cells for further analysis

mito.genes <- grep(pattern = "^MT-", x = rownames(x = EV@data), value = TRUE)
percent.mito <- colSums(EV@raw.data[mito.genes, ])/colSums(EV@raw.data)


EV <- AddMetaData(object = EV, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = EV, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = EV, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EV, gene1 = "nUMI", gene2 = "nGene")



EV <- FilterCells(object = EV, subset.names = c("nGene", "percent.mito"), 
                             low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.1))


# normalize
EV <- NormalizeData(object = EV, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# find variable gene
EV <- FindVariableGenes(object = EV, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
EV <- ScaleData(object = EV, vars.to.regress = c("nUMI", "percent.mito"))

########################################
# finding EV_mVenus positive cells

# ras+ by freebayes

EV_Ras<-data.frame(rownames=c("AACGTTGGTCAACATC","AAGACCTGTGCCTTGG","ACGTCAAGTGACAAAT","ACTTGTTAGATCTGCT","ATTGGTGCATGTCCTC","CAACCTCCAAACGTGG","CAGATCAGTAAATGAC","CCATTCGGTGCCTGTG","CCTAAAGTCGCGCCAA","CGACTTCGTATGAATG","CGCTGGACAGGTGGAT","CGTCACTGTACCGTAT","CTACGTCGTCTGCAAT","CTCCTAGAGGATGGTC","CTCGTCAAGCGTGAAC","CTGATAGCACTAAGTC","CTGTTTAAGATAGGAG","CTTCTCTTCTATGTGG","CTTGGCTCAGCCAATT","GAAACTCAGATATGCA","GATGCTATCCCTAACC","GCCTCTAAGCACCGCT","GCTTCCATCACAACGT","TGCTACCCATGTCCTC"
),Con1=rep("EV_RIS",24),Con2=rep("RIS",24)) # identified using step2_ras_snp.sh
rownames(EV_Ras)<-EV_Ras$rownames
EV_Ras<-EV_Ras[,-1]

EV_mVenus<-data.frame(t(data.frame(as.matrix(EV@data[grep("^mVenus$|^pLPC_puro$",rownames(EV@data)),]))))
EV_mVenus_pos<-data.frame(EV_mVenus[which(EV_mVenus$mVenus>0|EV_mVenus$pLPC_puro > 0),])
# finding EV_neo positive cells
EV_neo<-data.frame(EV@data[grep("^Neomycin$",rownames(EV@data)),])
EV_neo<-data.frame(EV_neo[which(EV_neo$EV.data.grep...Neomycin....rownames.EV.data....>0),,drop=F])
EV_neo_pos<-EV_neo
EV_Ras_pos<-EV_Ras
Reduce(intersect,list(rownames(EV_mVenus_pos),rownames(EV_neo),rownames(EV_Ras)))
intersect(rownames(EV_mVenus_pos),rownames(EV_neo))
intersect(rownames(EV_mVenus_pos),rownames(EV_Ras))
intersect(rownames(EV_neo),rownames(EV_Ras))

# mVenus only cell
EV_mVenus<-subset(EV_mVenus_pos, !(rownames(EV_mVenus_pos) %in% rownames(EV_neo)))
EV_mVenus<-subset(EV_mVenus, !(rownames(EV_mVenus) %in% rownames(EV_Ras)))

# of rasv12 only cell
EV_Ras<-subset(EV_Ras, !(rownames(EV_Ras) %in% rownames(EV_mVenus_pos)))
EV_Ras<-subset(EV_Ras, !(rownames(EV_Ras) %in% rownames(EV_neo_pos)))

EV_neo<-subset(EV_neo, !(rownames(EV_neo) %in% rownames(EV_mVenus_pos)))
# of neo only cells
EV_neo<-subset(EV_neo, !(rownames(EV_neo) %in% rownames(EV_Ras_pos)))

# neo and rasv12 only cells
EV_neo_ras<-subset(EV_neo_pos,rownames(EV_neo_pos) %in% rownames(EV_Ras_pos))
EV_neo_ras<-subset(EV_neo_ras,!(rownames(EV_neo_ras) %in% rownames(EV_mVenus_pos)))


EV_mVenus$Con1<-rep("EV_mVenus",dim(EV_mVenus)[1])
EV_neo$Con1<-rep("EV_RIS",dim(EV_neo)[1])
EV_mVenus$Con2<-rep("mVenus",dim(EV_mVenus)[1])
EV_neo$Con2<-rep("RIS",dim(EV_neo)[1])
EV_neo_ras$Con1<-rep("EV_RIS",dim(EV_neo_ras)[1])
EV_neo_ras$Con2<-rep("RIS",dim(EV_neo_ras)[1])

EV_neo_ras<-EV_neo_ras[,-1]
EV_mVenus<-EV_mVenus[,-c(1,2)]
EV_neo<-EV_neo[,-1]
EV_RIS_mVenus<-rbind(EV_mVenus,EV_neo,EV_Ras,EV_neo_ras)
#######################################################################################################

####################
# find cell labels that are mVenus+ and neo+

#######
EV_RIS_mVenus$Group<-paste(rep("EV_",dim(EV_RIS_mVenus)[1]),rownames(EV_RIS_mVenus),sep="")
dnMAML1_RIS_mVenus$Group<-paste(rep("dnMAML1_",dim(dnMAML1_RIS_mVenus)[1]),rownames(dnMAML1_RIS_mVenus),sep="")


passed<-data.frame(ID_Group=data.frame(rbind(EV_RIS_mVenus,dnMAML1_RIS_mVenus)$Group,celltype=data.frame(rbind(EV_RIS_mVenus,dnMAML1_RIS_mVenus))$Con1))
passed$Con2<-gsub("_.*","",passed$ID_Group.celltype)
rownames(passed)<-passed$ID_Group.rbind.EV_RIS_mVenus..dnMAML1_RIS_mVenus..Group
passed<-passed[,-1]




