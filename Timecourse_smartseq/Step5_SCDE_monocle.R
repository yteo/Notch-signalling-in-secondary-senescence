
library("scde")
library("ggplot2")

# IMR 90 time-course smart seq 
##SCDE
dat<-read.table("Compiled_downsampled_chr.txt")

cd<-dat
cd2<-data.frame(cd)
for (column in 1:224) {
  cd2[, column] <- as.numeric(as.character(cd2[, column]))
}

# omit genes that are never detected
cd2 <- cd2[rowSums(cd2)>0, ]

cd<-cd2[substr(row.names(cd2),1,4)=="ENSG",]
# omit cells with very poor coverage
ERCC<-cd2[substr(row.names(cd2),1,4)=="ERCC",]

cd <- cd[, colSums(cd)>5000]
cd<-apply(cd,2,function(cd) {storage.mode(cd) <- 'integer'; cd})

# omit genes that are not expressed in sufficient number of cells (more than 50)
cd <- cd[rowSums(cd>0)>50,]

cd_name<-merge(gene_name,cd,by.x="gene_id",by.y="row.names")
rownames(cd_name)<-cd_name$gene_id
cd_name<-cd_name[,-1]
write.table(cd_name,file="filtered_raw_IMR90.txt",quote=F,sep="\t")

con_filter <- factor(gsub("(Growing|day4|Sen|day2).*", "\\1", colnames(cd)), levels = c("Growing","day2","day4", "Sen"))
names(con_filter) <- colnames(cd)  
table(con_filter)


# calculate models
o.ifm <- scde.error.models(counts = cd, groups = con_filter, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

save(o.ifm,file="o.ifm_IMR90.RData")
#load("o.ifm_IMR90.RData")

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)




#Testing for differential expression IMR90


# define two groups of cells
groups <- factor(gsub("(Sen|Growing).*", "\\1", rownames(o.ifm)), levels  =  c("Sen", "Growing"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
#load("ediff.RData")
# convert z score to p value

p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
significant.genes <- which(p.values.adj<0.05)
length(significant.genes)


ord <- order(p.values.adj[significant.genes]) # order by p-value
de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
colnames(de) <- c("Lower bound","log2 fold change","Upper bound","p-value")

gene_name<-read.table("AllRNA_ID_name",header=T)
gene_name$gene_id<-sub("*\\..*", "",gene_name$gene_id)
de3<-cbind(de,rownames(de))
de4<-merge(de3,gene_name,by.x=c("rownames(de)"),by.y=c("gene_id"))


write.csv(de4,file="SCDE_DE2_newscde.csv",quote=FALSE)
################################################################################
# convert o.fpm to fpkm for monocle
# get expression magntiude estimates
o.fpm <- exp(scde.expression.magnitude(o.ifm, counts = cd))
gene_length<-read.table("exon_sizes",header=T)
rownames(gene_length)<-sub("*\\..*", "", rownames(gene_length))


fpkm<-o.fpm/gene_length[which(rownames(gene_length)%in%rownames(o.fpm)),]$Gene_length*1000

write.csv(fpkm,file="fpkm.csv",quote=F)
##################################################################################
# monocle
##########################################

# Loading data into monocle after SCDE
# Create a new CellDataSet

library(monocle)
fpkm_matrix <- read.csv("fpkm.csv",row.names=1)
gene_name<-read.table("AllRNA_ID_name",header=T)
gene_name$gene_id<-sub('\\..*', '', gene_name$gene_id)
gene_annotations<-data.frame(gene_id=rownames(fpkm_matrix))
gene_annotations<-merge(gene_annotations,gene_name,by="gene_id")
rownames(gene_annotations)<-gene_annotations$gene_id
gene_annotations$gene_id<-NULL
sample_sheet<-data.frame(cell_grp=sub("Sen1","Sen",sub("Sen3","Sen",sub("Sen2","Sen",sub("Growing2","Growing",sub("day2","Day2",sub("day4","Day4",gsub("\\_.*","",colnames(fpkm_matrix)))))))))
rownames(sample_sheet)<-colnames(fpkm_matrix)

# batch effect
sample_sheet2<-sample_sheet
sample_sheet2$batch<-rep("Batch1",dim(sample_sheet2)[1])
sample_sheet2$batch[grep("^day2",rownames(sample_sheet2))]<-"Batch2"
sample_sheet2$batch[grep(".*L003.*",rownames(sample_sheet2))]<-"Batch2"
sample_sheet2$batch[grep("^Sen3",rownames(sample_sheet2))]<-"Batch3"
sample_sheet2$pcolor[sample_sheet2$cell_grp=="Day2"] <- "palegreen4"
sample_sheet2$pcolor[sample_sheet2$cell_grp=="Day4"] <- "plum4"
sample_sheet2$pcolor[sample_sheet2$cell_grp=="Growing"] <- "skyblue"
sample_sheet2$pcolor[sample_sheet2$cell_grp=="Sen"] <- "red"

#fpkm_matrix[is.na(fpkm_matrix)] <- 0
pd <- new("AnnotatedDataFrame", data = sample_sheet2)
fd <- new("AnnotatedDataFrame",data=gene_annotations)
Timecourse <- newCellDataSet(as.matrix(fpkm_matrix), phenoData = pd, featureData = fd)
Timecourse<-estimateSizeFactors(Timecourse)
# choosing genes for filtering
Timecourse <- setOrderingFilter(Timecourse,rownames(de))
# reduce dimentionality
Timecourse<-reduceDimension(Timecourse, max_components=3)
# order cells
Timecourse <- orderCells(Timecourse, reverse=T)


plot_cell_trajectory(Timecourse, 1, 2, color="cell_grp")+geom_point(aes_string(color= "cell_grp", shape = "batch",),size=4)

plot_cell_trajectory(Timecourse,1,2,color="cell_grp",show_cell_names = F,cell_name_size=1)+geom_point(aes_string(color= "cell_grp"))+scale_color_manual(values=c("palegreen4","plum4","skyblue","red"))+ guides(colour = guide_legend(override.aes = list(size=5)))


save(Timecourse,file="Timecourse_monocle2.Rdata")


########################################
splitcells<-read.csv("split.csv",header=T)

# expression level of genes in monocle
###########################################################
# example to draw boxplots from fpkm--depending on data, here, i found the branching of sen into sentop and senbottom through monocle analysis. Hence, the following analysis consist of the split of fpkm table into the branch
box4<-function(id,genename)
{
  id<-deparse(substitute(id))
  gene<-deparse(substitute(genename))
  genename<-deparse(substitute(genename))

SenBottom_fpkm<-fpkm_matrix[, colnames(fpkm_matrix) %in% splitcells$Sen_17cells_individualbranch]
colnames(SenBottom_fpkm) <- paste("SenBottom", colnames(SenBottom_fpkm), sep = "_")

Senday24_fpkm<-fpkm_matrix[, colnames(fpkm_matrix) %in% splitcells$Sen_22cells_day24branch]
colnames(Senday24_fpkm) <- paste("SenTop", colnames(Senday24_fpkm), sep = "_")

Growing_split_fpkm<-fpkm_matrix[, grep("^Growing", colnames(fpkm_matrix))]
Sen_split_fpkm<-cbind(SenBottom_fpkm,Senday24_fpkm)

Day2_4_fpkm<-fpkm_matrix[, grep("^day2|day4", colnames(fpkm_matrix))]
fpkm_split<-cbind(Growing_split_fpkm,Sen_split_fpkm,Day2_4_fpkm)
fpkm_split<-log10(fpkm_split+1)
row<-rownames(fpkm_split)
fpkm_split<-as.data.frame(sapply(fpkm_split, as.numeric))
rownames(fpkm_split)<-row

library("ggplot2")
library("reshape")
gene.name<-id
# get expression magntiude estimates

toplot<-fpkm_split[row.names(fpkm_split)==gene.name,]
toplot<-melt(toplot)

toplot$con <- factor(gsub("(SenBottom|SenTop|Growing|day2|day4).*", "\\1", toplot$variable), levels = c("Growing","day2","day4", "SenTop","SenBottom"))

toplot$original<-toplot$variable

toplot<-toplot[which(is.finite(toplot$value)),]
#toplot$variable <- factor(toplot$variable , as.character(toplot$variable ))
ggplot(toplot,aes(x=con,y=value,fill=con))+geom_boxplot(outlier.shape=NA)+theme_bw()+ylab("log10 (FPKM+1)")+ggtitle(paste(genename))+theme( panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(text = element_text(size=12))+ theme(legend.position="bottom")+xlab("")+
    scale_fill_manual(values=c("skyblue","palegreen4","plum4","indianred","maroon"))+
    theme(legend.position="none")+scale_y_continuous(limits = c(0,4))

}
####################################

# example boxplot

box4(ENSG00000187498,COL4A1)











###################################
# SCDE of senescent branch

SenBottom<-cd[, colnames(cd) %in% splitcells$Sen_17cells_individualbranch]
colnames(SenBottom) <- paste("SenBottom", colnames(SenBottom), sep = "_")

Senday24<-cd[, colnames(cd) %in% splitcells$Sen_22cells_day24branch]
colnames(Senday24) <- paste("SenTop", colnames(Senday24), sep = "_")

Growing_split<-cbind(Growingtop,Growingbottom)
Sen_split<-cbind(SenBottom,Senday24)

#############################################################################################
######### sen SCDE split

con_sen <- factor(gsub("(SenBottom|SenTop).*", "\\1", colnames(Sen_split)))
names(con_sen) <- colnames(sen_split)  
table(con_sen)
o.ifm_sensplit <- scde.error.models(counts = Sen_split, groups = con_sen, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)


# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells_sen <- o.ifm_sensplit$corr.a > 0
table(valid.cells_sen)
o.ifm_sensplit <- o.ifm_sensplit[valid.cells_sen, ]

# estimate gene expression prior
o.prior_sen <- scde.expression.prior(models = o.ifm_sensplit, counts = Sen_split, length.out = 400, show.plot = FALSE)




#Testing for differential expression sen split


# define two groups_sen of cells
groups_sen <- factor(gsub("(SenBottom|SenTop).*", "\\1", rownames(o.ifm_sensplit)))
names(groups_sen) <- row.names(o.ifm_sensplit)
# run differential expression tests on all genes.
ediff_sen <- scde.expression.difference(o.ifm_sensplit, Sen_split, o.prior_sen, groups  =  groups_sen, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert z score to p value


p.values_sen <- 2*pnorm(abs(ediff_sen$Z),lower.tail=F) # 2-tailed p-value
p.values_sen.adj <- 2*pnorm(abs(ediff_sen$cZ),lower.tail=F) # Adjusted to control for FDR
significant.genes_sen <- which(p.values_sen.adj<0.05)
length(significant.genes_sen)
##
# GSEA
allgene<-which(p.values_sen.adj<1.1)
ord_all<-order(p.values_sen.adj[allgene])
sen_split_GSEA <- data.frame(cbind(ediff_sen[allgene,1:3],p.values_sen.adj[allgene])[ord_all,])
colnames(sen_split_GSEA) <- c("Lower bound","log2 fold change","Upper bound","p-value")
gene_name<-read.table("AllRNA_ID_name",header=T)
gene_name$gene_id<-sub("*\\..*", "",gene_name$gene_id)
sen_split_GSEA<-merge(gene_name,sen_split_GSEA,by.y="row.names",by.x=c("gene_id"))

write.csv(sen_split_GSEA,file="sen_split_GSEA.csv",quote=F)
##
ord <- order(p.values_sen.adj[significant.genes_sen]) # order by p-value
de_sen <- cbind(ediff_sen[significant.genes_sen,1:3],p.values_sen.adj[significant.genes_sen])[ord,]
colnames(de_sen) <- c("Lower bound","log2 fold change","Upper bound","p-value")

gene_name<-read.table("AllRNA_ID_name",header=T)
gene_name$gene_id<-sub("*\\..*", "",gene_name$gene_id)
de3_sen<-cbind(de_sen,rownames(de_sen))
de4_sen<-merge(de3_sen,gene_name,by.x=c("rownames(de_sen)"),by.y=c("gene_id"))
