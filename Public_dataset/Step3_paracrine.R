library(limma)
#RMA.txt obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41318
dat<-read.table("RMA.txt",header=T,row.names=1)
design <- cbind(Control=1,Paracrine=c(0,0,0,1,1,1))
fit<-lmFit(dat,design)
fit <- eBayes(fit)
FC<-topTable(fit,coef=2,number=Inf)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)


library(annotate)
library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"
ID     <- rownames(FC)
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))


FC$Symbol<-Symbol

FC<-FC[which(FC$Symbol!="NA"),]
gsea<-data.frame(FC$Symbol,FC$logFC)
write.table(gsea,file="Paracrine_Growing.rnk",sep="\t",quote=F,row.names = F)


FC_sig<-FC[FC$adj.P.Val<0.05,]

write.table(FC_sig,file="Paracrine_ControlGrowing.txt",sep="\t")


##############################

# OIS vs paracrine

library(limma)
dat<-read.table("OIS_paracrine_RMA.txt",header=T,row.names=1)
design <- cbind(OIS=1,Paracrine=c(0,0,0,1,1,1))
fit<-lmFit(dat,design)
fit <- eBayes(fit)
FC<-topTable(fit,coef=2,number=Inf)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

convert<-read.table("OIS_paracrine_geneID2.txt")
FC2<-merge(FC,convert,by.x="row.names",by.y="V2")
library(annotate)
library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"
ID     <- rownames(FC)
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))


FC$Symbol<-Symbol

FC<-FC[which(FC$Symbol!="NA"),]

FC_sig<-FC[FC$adj.P.Val<0.05,]
FC_sig2<-FC2[FC2$adj.P.Val<0.05,]
write.table(FC_sig,file="Paracrine_OIS.txt",sep="\t",quote=F)
write.table(FC_sig2,file="Paracrine_OIS_correctID.txt",sep="\t",quote=F)
gsea<-data.frame(FC$Symbol,FC$logFC)
write.table(gsea,file="Paracrine_OIS.rnk",sep="\t",quote=F,row.names = F)

