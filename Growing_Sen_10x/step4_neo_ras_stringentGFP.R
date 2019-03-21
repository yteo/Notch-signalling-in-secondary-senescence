library(Seurat)
library(dplyr)
library(Matrix)

# Load the imr90sen1 dataset
imr90sen1.data <- Read10X(data.dir ="IMR90_Senrep1_merged_2_0_1\\outs\\filtered_gene_bc_matrices\\hg19_neo_GFP_puro\\")
#x = as.matrix(x = imr90sen1.data)

imr90sen1 <- CreateSeuratObject(raw.data = imr90sen1.data,min.cells=3,min.genes=200,project="10X_IMR90sen1")


#####################
# basic QC and selecting cells for further analysis

mito.genes <- grep(pattern = "^MT-", x = rownames(x = imr90sen1@data), value = TRUE)
percent.mito <- colSums(imr90sen1@raw.data[mito.genes, ])/colSums(imr90sen1@raw.data)

imr90sen1 <- AddMetaData(object = imr90sen1, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = imr90sen1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = imr90sen1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = imr90sen1, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts less than 2500 and MT greater than 0.15
imr90sen1 <- FilterCells(object = imr90sen1, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.15))


# normalize
imr90sen1 <- NormalizeData(object = imr90sen1, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

# find variable gene
imr90sen1 <- FindVariableGenes(object = imr90sen1, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
imr90sen1 <- ScaleData(object = imr90sen1, vars.to.regress = c("nUMI", "percent.mito"))

########################################
# finding sen1_GFP positive cells
sen1_Ras<-data.frame(rownames=c("AAGGTTCAGACAGGCT","ACACCGGTCTCAAACG","ACACTGAGTTCCACAA","ACACTGATCCAAATGC","ACATACGTCATTATCC","ACATCAGTCGCATGGC","ACCCACTAGACTAGAT","ACGGAGAAGAGGGATA","ACGGAGACACACATGT","ACGGGCTAGCCAGTTT","ACTTTCAGTTAGATGA","AGCCTAAGTCGCCATG","AGCTTGACAGTAACGG","AGGCCACAGGCTACGA","AGGGTGACAGCGATCC","ATCTGCCGTGTCAATC","ATTATCCTCCTAAGTG","ATTATCCTCTGTGCAA","CAAGGCCTCCCATTTA","CAAGTTGTCTCCGGTT","CACACTCTCTTGGGTA","CACATAGTCTGCAAGT","CAGCCGAGTGCAACTT","CAGCCGATCCTATGTT","CAGGTGCTCCTGCTTG","CAGTAACAGCTACCTA","CATCAAGGTCGAGATG","CCACTACTCAGCGACC","CCTACACCAAGAGGCT","CCTATTAGTCGGCTCA","CCTCTGAAGGGAGTAA","CCTTCGAGTCTCGTTC","CGCCAAGTCTAACTTC","CGGACTGGTAGAGTGC","CGGAGTCAGGCCGAAT","CGTAGCGGTAATCGTC","CGTCTACCACCACCAG","CTCAGAACAAGACGTG","CTCTAATCACTTAACG","CTCTGGTGTCAACTGT","CTGAAACTCTACTATC","CTGATCCGTGCCTGCA","CTTCTCTTCAGTTGAC","GAACATCGTCTAGGTT","GAACCTAGTTTGCATG","GAACGGAAGTGTACGG","GACAGAGGTTGTACAC","GAGCAGAAGCTTATCG","GAGTCCGTCTTGCATT","GATTCAGCACATTCGA","GCATGATAGTCTCAAC","GCGACCATCCCAAGTA","GCGGGTTGTGCAACTT","GCTTGAAAGCTTCGCG","GGATGTTCAACTGCTA","GGCTGGTGTGGTACAG","GGTGAAGCAATCGGTT","GGTGTTATCACAAACC","GTAGGCCCAACACGCC","GTAGGCCTCCTGCTTG","GTCTTCGAGAAACGAG","GTTCTCGCACTCGACG","TACACGACAAGGACTG","TCAGATGAGTCTTGCA","TCCACACCACCAACCG","TCGCGAGAGCTGATAA","TCGCGAGTCGGACAAG","TCGTACCGTAATAGCA","TCGTAGAGTTACAGAA","TGCGGGTGTACCAGTT","TGTATTCGTTGTTTGG","TGTATTCTCGCAGGCT","TTCCCAGCAAGCCGTC"
),Con1=rep("sen1_RIS",73),Con2=rep("RIS",73)) # identified using step3_ras_snp.sh
rownames(sen1_Ras)<-sen1_Ras$rownames
sen1_Ras<-sen1_Ras[,-1]
sen1_GFP<-data.frame(t(data.frame(as.matrix(imr90sen1@data[grep("^GFP$|^puro$",rownames(imr90sen1@data)),]))))
sen1_GFP_pos<-data.frame(sen1_GFP[which(sen1_GFP$GFP>0.3|sen1_GFP$puro > 0.3),])
# finding sen1_neo positive cells
sen1_neo<-data.frame(imr90sen1@data[grep("^ras_neo$",rownames(imr90sen1@data)),])
sen1_neo<-data.frame(sen1_neo[which(sen1_neo$imr90sen1.data.grep...ras_neo....rownames.imr90sen1.data.....>0),,drop=F])
sen1_neo_pos<-sen1_neo
sen1_Ras_pos<-sen1_Ras
Reduce(intersect,list(rownames(sen1_GFP_pos),rownames(sen1_neo),rownames(sen1_Ras)))
intersect(rownames(sen1_GFP_pos),rownames(sen1_neo))
intersect(rownames(sen1_GFP_pos),rownames(sen1_Ras))
intersect(rownames(sen1_neo),rownames(sen1_Ras))

# gfp only cell
sen1_GFP<-subset(sen1_GFP_pos, !(rownames(sen1_GFP_pos) %in% rownames(sen1_neo)))
sen1_GFP<-subset(sen1_GFP, !(rownames(sen1_GFP) %in% rownames(sen1_Ras)))

# of rasv12 only cell
sen1_Ras<-subset(sen1_Ras, !(rownames(sen1_Ras) %in% rownames(sen1_GFP_pos)))
sen1_Ras<-subset(sen1_Ras, !(rownames(sen1_Ras) %in% rownames(sen1_neo_pos)))

sen1_neo<-subset(sen1_neo, !(rownames(sen1_neo) %in% rownames(sen1_GFP_pos)))
# of neo only cells
sen1_neo<-subset(sen1_neo, !(rownames(sen1_neo) %in% rownames(sen1_Ras_pos)))

# neo and rasv12 only cells
sen1_neo_ras<-subset(sen1_neo_pos,rownames(sen1_neo_pos) %in% rownames(sen1_Ras_pos))
sen1_neo_ras<-subset(sen1_neo_ras,!(rownames(sen1_neo_ras) %in% rownames(sen1_GFP_pos)))


sen1_GFP$Con1<-rep("sen1_GFP",dim(sen1_GFP)[1])
sen1_neo$Con1<-rep("sen1_RIS",dim(sen1_neo)[1])
sen1_GFP$Con2<-rep("GFP",dim(sen1_GFP)[1])
sen1_neo$Con2<-rep("RIS",dim(sen1_neo)[1])
sen1_neo_ras$Con1<-rep("sen1_RIS",dim(sen1_neo_ras)[1])
sen1_neo_ras$Con2<-rep("RIS",dim(sen1_neo_ras)[1])

sen1_neo_ras<-sen1_neo_ras[,-1]
sen1_GFP<-sen1_GFP[,-c(1,2)]
sen1_neo<-sen1_neo[,-1]
sen1_RIS_GFP<-rbind(sen1_GFP,sen1_neo,sen1_Ras,sen1_neo_ras)
########################################################################################################################################################
# growing 1
library(Seurat)
library(dplyr)
library(Matrix)
# Load the imr90growing1 dataset
imr90growing1.data <- Read10X(data.dir ="IMR90_growingrep1_merged_2_0_1\\outs\\filtered_gene_bc_matrices\\hg19_neo_GFP_puro\\")


imr90growing1 <- CreateSeuratObject(raw.data = imr90growing1.data,min.cells=3,min.genes=200,project="10X_IMR90growing1")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = imr90growing1@data), value = TRUE)
percent.mito <- colSums(imr90growing1@raw.data[mito.genes, ])/colSums(imr90growing1@raw.data)


imr90growing1 <- AddMetaData(object = imr90growing1, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = imr90growing1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = imr90growing1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = imr90growing1, gene1 = "nUMI", gene2 = "nGene")


imr90growing1 <- FilterCells(object = imr90growing1, subset.names = c("nGene", "percent.mito"), 
                             low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.15))


# normalize
imr90growing1 <- NormalizeData(object = imr90growing1, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# find variable gene
imr90growing1 <- FindVariableGenes(object = imr90growing1, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
imr90growing1 <- ScaleData(object = imr90growing1, vars.to.regress = c("nUMI", "percent.mito"))

########################################
# finding growing1_GFP positive cells

# ras+ by freebayes

growing1_Ras<-data.frame(rownames=c("AGAGTGGAGTCCCACG","AGTGGGAGTTGTCGCG","CAAGTTGAGTTACGGG","GGCTCGAGTCACACGC","TACGGATTCATCATTC","TACTTGTAGATGGCGT","TCATTACAGGGTCGAT","TGACAACTCATTTGGG","TTGTAGGAGCCAGTTT"
),Con1=rep("growing1_RIS",9),Con2=rep("RIS",9))
rownames(growing1_Ras)<-growing1_Ras$rownames
growing1_Ras<-growing1_Ras[,-1]

growing1_GFP<-data.frame(t(data.frame(as.matrix(imr90growing1@data[grep("^GFP$|^puro$",rownames(imr90growing1@data)),]))))
growing1_GFP_pos<-data.frame(growing1_GFP[which(growing1_GFP$GFP>0.3|growing1_GFP$puro > 0.3),])
# finding growing1_neo positive cells
growing1_neo<-data.frame(imr90growing1@data[grep("^ras_neo$",rownames(imr90growing1@data)),])
growing1_neo<-data.frame(growing1_neo[which(growing1_neo$imr90growing1.data.grep...ras_neo....rownames.imr90growing1.data....>0),,drop=F])
growing1_neo_pos<-growing1_neo
growing1_Ras_pos<-growing1_Ras
Reduce(intersect,list(rownames(growing1_GFP_pos),rownames(growing1_neo),rownames(growing1_Ras)))
intersect(rownames(growing1_GFP_pos),rownames(growing1_neo))
intersect(rownames(growing1_GFP_pos),rownames(growing1_Ras))
intersect(rownames(growing1_neo),rownames(growing1_Ras))

# gfp only cell
growing1_GFP<-subset(growing1_GFP_pos, !(rownames(growing1_GFP_pos) %in% rownames(growing1_neo)))
growing1_GFP<-subset(growing1_GFP, !(rownames(growing1_GFP) %in% rownames(growing1_Ras)))

# of rasv12 only cell
growing1_Ras<-subset(growing1_Ras, !(rownames(growing1_Ras) %in% rownames(growing1_GFP_pos)))
growing1_Ras<-subset(growing1_Ras, !(rownames(growing1_Ras) %in% rownames(growing1_neo_pos)))

growing1_neo<-subset(growing1_neo, !(rownames(growing1_neo) %in% rownames(growing1_GFP_pos)))
# of neo only cells
growing1_neo<-subset(growing1_neo, !(rownames(growing1_neo) %in% rownames(growing1_Ras_pos)))

# neo and rasv12 only cells
growing1_neo_ras<-subset(growing1_neo_pos,rownames(growing1_neo_pos) %in% rownames(growing1_Ras_pos))
growing1_neo_ras<-subset(growing1_neo_ras,!(rownames(growing1_neo_ras) %in% rownames(growing1_GFP_pos)))


growing1_GFP$Con1<-rep("growing1_GFP",dim(growing1_GFP)[1])
growing1_neo$Con1<-rep("growing1_RIS",dim(growing1_neo)[1])
growing1_GFP$Con2<-rep("GFP",dim(growing1_GFP)[1])
growing1_neo$Con2<-rep("RIS",dim(growing1_neo)[1])
growing1_neo_ras$Con1<-rep("growing1_RIS",dim(growing1_neo_ras)[1])
growing1_neo_ras$Con2<-rep("RIS",dim(growing1_neo_ras)[1])

growing1_neo_ras<-growing1_neo_ras[,-1]
growing1_GFP<-growing1_GFP[,-c(1,2)]
growing1_neo<-growing1_neo[,-1]
growing1_RIS_GFP<-rbind(growing1_GFP,growing1_neo,growing1_Ras,growing1_neo_ras)
#######################################################################################################

# growing 2
library(Seurat)
library(dplyr)
library(Matrix)
# Load the imr90growing2 dataset
imr90growing2.data <- Read10X(data.dir ="IMR90_growingrep2_merged_2_0_1\\outs\\filtered_gene_bc_matrices\\hg19_neo_GFP_puro\\")
#x = as.matrix(x = imr90growing2.data)

imr90growing2 <- CreateSeuratObject(raw.data = imr90growing2.data,min.cells=3,min.genes=200,project="10X_IMR90growing2")


#####################
# basic QC and selecting cells for further analysis

mito.genes <- grep(pattern = "^MT-", x = rownames(x = imr90growing2@data), value = TRUE)
percent.mito <- colSums(imr90growing2@raw.data[mito.genes, ])/colSums(imr90growing2@raw.data)

imr90growing2 <- AddMetaData(object = imr90growing2, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = imr90growing2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = imr90growing2, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = imr90growing2, gene1 = "nUMI", gene2 = "nGene")

imr90growing2 <- FilterCells(object = imr90growing2, subset.names = c("nGene", "percent.mito"), 
                             low.thresholds = c(2500, -Inf), high.thresholds = c(Inf, 0.15))


# normalize
imr90growing2 <- NormalizeData(object = imr90growing2, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# find variable gene
imr90growing2 <- FindVariableGenes(object = imr90growing2, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scaling and regress out effects
imr90growing2 <- ScaleData(object = imr90growing2, vars.to.regress = c("nUMI", "percent.mito"))

########################################
# finding growing2_GFP positive cells

# ras+ by freebayes


growing2_Ras<-data.frame(rownames=c("AAGGAGCGTTGATTGC","ATTGGTGTCTTACCGC","CAGGTGCAGGATGGTC","CCAATCCTCAGAGGTG","CCCAGTTCAAGCGTAG","GAGGTGAGTGCGCTTG","GCTGCTTAGGGCATGT","GTATCTTAGACTGTAA","TGTGGTAAGAGGTTAT","TTCTCAATCAACGGGA","TTTGTCACAATCTGCA"
),Con1=rep("growing2_RIS",11),Con2=rep("RIS",11))
rownames(growing2_Ras)<-growing2_Ras$rownames
growing2_Ras<-growing2_Ras[,-1]
growing2_GFP<-data.frame(t(data.frame(as.matrix(imr90growing2@data[grep("^GFP$|^puro$",rownames(imr90growing2@data)),]))))
growing2_GFP_pos<-data.frame(growing2_GFP[which(growing2_GFP$GFP>0.3|growing2_GFP$puro > 0.3),])
# finding growing2_neo positive cells
growing2_neo<-data.frame(imr90growing2@data[grep("^ras_neo$",rownames(imr90growing2@data)),])
growing2_neo<-data.frame(growing2_neo[which(growing2_neo$imr90growing2.data.grep...ras_neo....rownames.imr90growing2.data....>0),,drop=F])
growing2_neo_pos<-growing2_neo
growing2_Ras_pos<-growing2_Ras
Reduce(intersect,list(rownames(growing2_GFP_pos),rownames(growing2_neo),rownames(growing2_Ras)))
intersect(rownames(growing2_GFP_pos),rownames(growing2_neo))
intersect(rownames(growing2_GFP_pos),rownames(growing2_Ras))
intersect(rownames(growing2_neo),rownames(growing2_Ras))

# gfp only cell
growing2_GFP<-subset(growing2_GFP_pos, !(rownames(growing2_GFP_pos) %in% rownames(growing2_neo)))
growing2_GFP<-subset(growing2_GFP, !(rownames(growing2_GFP) %in% rownames(growing2_Ras)))

# of rasv12 only cell
growing2_Ras<-subset(growing2_Ras, !(rownames(growing2_Ras) %in% rownames(growing2_GFP_pos)))
growing2_Ras<-subset(growing2_Ras, !(rownames(growing2_Ras) %in% rownames(growing2_neo_pos)))

growing2_neo<-subset(growing2_neo, !(rownames(growing2_neo) %in% rownames(growing2_GFP_pos)))
# of neo only cells
growing2_neo<-subset(growing2_neo, !(rownames(growing2_neo) %in% rownames(growing2_Ras_pos)))

# neo and rasv12 only cells
growing2_neo_ras<-subset(growing2_neo_pos,rownames(growing2_neo_pos) %in% rownames(growing2_Ras_pos))
growing2_neo_ras<-subset(growing2_neo_ras,!(rownames(growing2_neo_ras) %in% rownames(growing2_GFP_pos)))


growing2_GFP$Con1<-rep("growing2_GFP",dim(growing2_GFP)[1])
growing2_neo$Con1<-rep("growing2_RIS",dim(growing2_neo)[1])
growing2_GFP$Con2<-rep("GFP",dim(growing2_GFP)[1])
growing2_neo$Con2<-rep("RIS",dim(growing2_neo)[1])
growing2_neo_ras$Con1<-rep("growing2_RIS",dim(growing2_neo_ras)[1])
growing2_neo_ras$Con2<-rep("RIS",dim(growing2_neo_ras)[1])

growing2_neo_ras<-growing2_neo_ras[,-1]
growing2_GFP<-growing2_GFP[,-c(1,2)]
growing2_neo<-growing2_neo[,-1]
growing2_RIS_GFP<-rbind(growing2_GFP,growing2_neo,growing2_Ras,growing2_neo_ras)

####################
# find cell labels that are GFP+ and neo+

#######
growing1_RIS_GFP$Group<-paste(rownames(growing1_RIS_GFP),rep("-3",dim(growing1_RIS_GFP)[1]),sep="")
growing2_RIS_GFP$Group<-paste(rownames(growing2_RIS_GFP),rep("-4",dim(growing2_RIS_GFP)[1]),sep="")
sen1_RIS_GFP$Group<-paste(rownames(sen1_RIS_GFP),rep("-1",dim(sen1_RIS_GFP)[1]),sep="")
sen2_RIS_GFP$Group<-paste(rownames(sen2_RIS_GFP),rep("-2",dim(sen2_RIS_GFP)[1]),sep="")


passed<-data.frame(ID_Group=data.frame(rbind(growing1_RIS_GFP,growing2_RIS_GFP,sen1_RIS_GFP,sen2_RIS_GFP)$Group,celltype=data.frame(rbind(growing1_RIS_GFP,growing2_RIS_GFP,sen1_RIS_GFP,sen2_RIS_GFP))$Con1))
passed$Con2<-gsub("_.*","",passed$ID_Group.celltype)
rownames(passed)<-passed$ID_Group.rbind.growing1_RIS_GFP..growing2_RIS_GFP..sen1_RIS_GFP..sen2_RIS_GFP..Group
passed<-passed[,-1]
passed$Con3<-rep("Growing",dim(passed)[1])
passed$Con3[passed$ID_Group.celltype=="sen1_RIS"]<-"RIS"
passed$Con3[passed$ID_Group.celltype=="sen2_RIS"]<-"RIS"
passed$Con3[passed$ID_Group.celltype=="sen1_GFP"]<-"GFP"
passed$Con3[passed$ID_Group.celltype=="sen2_GFP"]<-"GFP"



write.table(passed,file="passedcell_GFP_RIS_stringentGFP.txt",sep="\t",quote=F,row.names = T,col.names = T)