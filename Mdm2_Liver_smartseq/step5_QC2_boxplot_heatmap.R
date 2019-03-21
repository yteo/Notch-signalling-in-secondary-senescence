library("ggplot2")

library(matrixStats)
library(tidyr)
dat<-read.table("Compiled_downsampled50000_chr.txt")

dat<-spread(dat,V1,V3)
rownames(dat)<-dat$V2
dat<-dat[,-1]


dat<-dat[apply(dat, 1, function(x) !all(x==0)),]
expr.count <- dat[grep("^ERCC",row.names(dat),invert=TRUE),]




total.sum.counts <- apply(dat,2,sum)


#Calculate the number of genes with at least 1 read in each sample (cell)


genes.with.reads <- sapply(1:ncol(expr.count),function(x){sum (!expr.count[,x] == 0)})
names(genes.with.reads) <- colnames(expr.count)



#Plot num genes with at least 1 read vs total gene count


con <- factor(gsub("(6_1|6_2).*", "\\1", colnames(dat)), levels = c("6_1","6_2"))

ggplot()+geom_point(aes(x=total.sum.counts,y=genes.with.reads,color=con))+xlab("Total gene count")+ylab("Num of genes with at least one read")+theme_bw()+scale_color_manual(name="",labels=c("Induced","Uninduced"),values=c("black","grey"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank(),plot.title = element_text(hjust = 1)) +ggtitle("64 cells")+geom_hline(yintercept=500,linetype="dashed")+
  geom_vline(xintercept=20000,linetype="dashed")

name<-as.data.frame(total.sum.counts)
name$con<-row.names(name)
passed<-name[( total.sum.counts > 20000 & genes.with.reads>500 ),]$con
dat_passed<-expr.count[,colnames(expr.count) %in% passed]



#####################


library(SC3)
library(tidyr)
library(scater)
library(scran)

dat<-dat_passed

colnames(dat)<-gsub("X","",colnames(dat))


cd_sc3<-dat[apply(dat, 1, function(x) !all(x==0)),]



dat2<-as.matrix(cd_sc3)
cluster<-data.frame(row.names=colnames(dat2),cluster=gsub("_.*","",colnames(dat2)))
sce <- SingleCellExperiment(assays=list(counts=dat2),colData=cluster)

is.mito <- grepl("^mt\\.", rownames(sce))
summary(is.mito)
sce <- calculateQCMetrics(sce,feature_controls=list( Mt=is.mito))

head(colnames(colData(sce)))

sce$plates <- gsub("(6_1|6_2).*","\\1",colnames(sce))
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")




libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE, batch=sce$plates)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE, batch=sce$plates)


sce <- sce[,!(libsize.drop | feature.drop )]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           Remaining=ncol(sce))


passed<-colnames(sce)
#write.table(passed,file="Notdownsampled_sce_QC2passed.txt",quote=F,col.names=F,row.names = F)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize

plotQC(sce, type = "highest-expression", n=50) + fontsize


ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

############

num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))



to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)

sce <- computeSumFactors(sce)
summary(sizeFactors(sce))


plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
     ylab="Library size (millions)", xlab="Size factor")





sce <- normalize(sce)

var.fit.nospike <- trendVar(sce, parametric=TRUE, use.spikes=FALSE, span=0.2)
var.out.nospike <- decomposeVar(sce, var.fit.nospike)

plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
     xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)


Denoising expression values using PCA

sce <- denoisePCA(sce, technical=var.fit.nospike$trend) 
dim(reducedDim(sce, "PCA")) 




exon2<-read.table("sce_mdm2Exon56mapped.lst")
mapstat<-read.table("mappingstat2.txt")

mdmneg<-exon2[exon2$V2==0,]$V1

mdmpos<-exon2[exon2$V2!=0,]$V1

# exprs returns log normalized count
exprs_mat <- exprs(sce)


mdmnegexpr<-exprs_mat[,which(colnames(exprs_mat) %in% mdmneg)]
mdmneg<-colnames(mdmnegexpr)
colnames(mdmnegexpr)<-paste("MdmNeg",colnames(mdmnegexpr),sep="_")

mdmposexpr<-exprs_mat[,which(colnames(exprs_mat) %in% mdmpos)]
mdmpos<-colnames(mdmposexpr)
colnames(mdmposexpr)<-paste("MdmPos",colnames(mdmposexpr),sep="_")

exprs_mat2<-cbind(mdmposexpr,mdmnegexpr)
'%!in%' <- function(x,y)!('%in%'(x,y))
exprs_mat3<-exprs_mat[,colnames(exprs_mat) %!in% gsub(".*6_","6_",colnames(exprs_mat2))]
colnames(exprs_mat3)<-paste("WT_",colnames(exprs_mat3),sep="")

all<-cbind(exprs_mat2,exprs_mat3)


ensembl<-read.table("modified_ENSM_genename_mm10.txt")
all<-merge(ensembl,all,by.x="V1",by.y="row.names")


colnames(all)[2]<-"name"

rownames(all)<-make.names(all$name,unique=T)
all<-all[,-1]
all<-all[,-1]



exon2_neg<-exon2[exon2$V1 %in% mdmneg,]
exon2_neg$exon<-"Mdm2-"
exon2_pos<-exon2[exon2$V1 %in% mdmpos,]
exon2_pos$exon<-"Mdm2+"
exon3<-rbind(exon2_neg,exon2_pos)
exon3$category[exon3$V2==0]<-"0"
exon3$category[(exon3$V2>0) & (exon3$V2<=500)]<-"1-500"
exon3$category[(exon3$V2>500) & (exon3$V2<=1000)]<-"501-1000"
exon3$category[(exon3$V2>1000) ]<-">1000"
exon3$category<-factor(exon3$category,levels = c("0","1-500","501-1000",">1000"))
ggplot()+geom_histogram(data=exon3,aes(x=category,fill="exon"),stat="count")+ylab("Number of cells")+xlab("Number of reads mapping to exon5/6 of Mdm2 gene")+theme_bw()

###############################


# boxplot and heatmap
# boxplot
####################
box<-function(gene){
  gene1<-t(data.frame(all[grep(paste("^",gene,"$",sep=""),rownames(all)),]))
  colnames(gene1)<-"exp"
  gene1<-data.frame(gene1)
  gene1$cluster<-gsub("_.*","",as.character(rownames(gene1)))
  
  
  ggplot(data=gene1,aes(x=cluster,y=exp))+geom_boxplot(outlier.colour=NA)+ylab("normalized expression)")+theme_bw()+ggtitle(paste(gene))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
box("Cdkn1a")

###################################


exprs<-exprs_mat2
colnames(exprs)<-factor(gsub("_.*","",colnames(exprs)),levels=c("MdmNeg","MdmPos"))
exprs<-exprs[,order(colnames(exprs))]
toMatch<-c("Maml1", "Rfng", "Dvl3", "Psenen", "Jag2", "Snw1", "Rbpj")
matches <-exprs[grep(paste("^",toMatch,"$",sep="",collapse="|"), 
                     rownames(exprs)),]
matches[matches >0] <- 1
library(gplots)
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(matches) , cluster_rows=FALSE, cluster_cols=FALSE, col=colfunc(15),dendrogram="none", density.info="none",trace="none",sepcolor="black",sepwidth=c(0.1,0.1))


colnames(matches)<-gsub("\\..*","",colnames(matches))
#number of cells that express transcript in Mdm2pos
sum(rowSums(matches[,grep("MdmPos",colnames(matches))]))
# number of cells that express transcript in Mdm2neg
sum(rowSums(matches[,grep("MdmNeg",colnames(matches))]))

#number of cells that do not express transcript in Mdm2pos
(25*7)-sum(rowSums(matches[,grep("MdmPos",colnames(matches))]))


(21*7)-sum(rowSums(matches[,grep("MdmNeg",colnames(matches))]))

tgf<-matrix(c(140,35,144,3),nrow=2,ncol=2)
colnames(tgf)<-c("Mdmpos","Mdmneg")
rownames(tgf)<-c("Expressed","Nonexpressed")
fisher.test(tgf)
#p-value = 1.447e-07


######## housekeeping genes

exprs<-exprs_mat2

colnames(exprs)<-factor(gsub("_.*","",colnames(exprs)),levels=c("MdmNeg","MdmPos"))
exprs<-exprs[,order(colnames(exprs))]

rownames(exprs)<-factor(rownames(exprs),levels=c("Tubb4b","Actb","Cdkn1a"))

toMatch<-c("Tubb4b","Actb","Cdkn1a")
matches <-exprs[grep(paste("^",toMatch,"$",sep="",collapse="|"), 
                     rownames(exprs)),]
library(gplots)
library(pheatmap)
colfunc <- colorRampPalette(c("white", "red"))

heatmap.2(matches , Rowv=FALSE, Colv=FALSE, col=colfunc(15),dendrogram="none", density.info="none",trace="none",sepcolor="black",sepwidth=c(0.1,0.1))

pheatmap(matches , cluster_rows=FALSE, cluster_cols=FALSE, col=colfunc(15),dendrogram="none", density.info="none",trace="none",sepcolor="black",sepwidth=c(0.1,0.1))


################################################################################################################
