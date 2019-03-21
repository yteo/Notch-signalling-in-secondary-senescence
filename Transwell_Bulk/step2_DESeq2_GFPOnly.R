library(DESeq2)

directory <- "./"
sampleFiles <- grep(".*_htseq.txt",list.files(directory),value=TRUE)
sampleCondition <- gsub("(_).*","",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

dds <-DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                 directory = directory,
                                 design= ~  condition)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]

dds <- DESeq(dds)



condition<-dds$condition



# Plot dispersions

plotDispEsts(dds, main="Dispersion plot")


# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))

library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)

heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")


rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
                                                                                      terldt = list(levels(fac)), rep = FALSE)))
}

rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))



########################################
# individual comparison
dds<-DESeq(dds)
res <- results(dds)
sizeFactors(dds)
res<-data.frame(res)

##############
# replace ENSEMBL ID with gene name
ID_name<-read.table("ENSEMBLID_name.txt",header=F,sep="\t")

ID_name$V1<-gsub("\\..*","",ID_name$V1)
rownames(res)<-gsub("\\..*","",rownames(res))
res<-merge(res,ID_name,by.x="row.names",by.y="V1")


# get counts
normCount<-log2(counts(dds, normalized=TRUE)+1)
rownames(normCount)<-gsub("\\..*","",rownames(normCount))
normCount<-merge(normCount,ID_name,by.x="row.names",by.y="V1")
colnames(normCount)<-gsub("(_).*","",colnames(normCount))


cnt<-counts(dds)

write.table(cnt,"Transwell_HTSeqCount.txt",sep="\t",quote=F)




res<-res[complete.cases(res), ]
write.table(data.frame(gene=res$V2,log2FC=res$log2FoldChange),file="GFPOut_In.rnk",sep="\t",quote=F,row.names = F)
write.table(res,file="GFPOut_In.txt",sep="\t",quote=F,row.names = F)


###############################################
# heatmap
features.plot<-c("^TGFBI$","^TGFB1$","^COL6A1$","^COL6A3$","^COL5A3$","^COL18A1$","^E2F7$")
matches <- normCount[grepl(paste(features.plot, collapse="|"), normCount$V2), ]

matches<-matches[,-1]
rownames(matches)<-matches$V2
matches<-matches[,-7]
library(RColorBrewer)
colnames(matches)<-gsub("\\..*","",colnames(matches))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(colnames(matches)))])


library(gplots)

heatmap.2(as.matrix(matches), key=T, trace="none",
          col=bluered(16),scale="row",
           ColSideColors=mycols[colnames(matches)],
          margin=c(10, 10), main="Gene expression")


################################################

# webgestalt
# ratio of enrichment
path=data.frame( source=rep("Wikipathway",8),pathway=c("Lung fibrosis","miRNA targets in ECM and membrane receptors","Matrix Metalloproteinases","PI3K-Akt Signaling Pathway","Differentiation Pathway","TGF-beta Receptor Signaling","Senescence and Autophagy in Cancer","Focal Adhesion-PI3K-Akt-mTOR-signaling pathway"),ratioenrichment=
                   c(9.89,10.38,12.71,2.85,8.9,5.52,3.75,2.64))

kegg=data.frame( source=rep("Kegg",7),pathway=c("Cytokine-cytokine receptor interaction","ECM-receptor interaction","Hematopoietic cell lineage","Proteoglycans in cancer","PI3K-Akt signaling pathway","MicroRNAs in cancer","Focal adhesion"),ratioenrichment=
                   c(6.22,6.08,8.41,4.56,2.86,3.32,3.13))


path<-data.frame(rbind(path,kegg))
library(ggplot2)
ggplot(data=path,aes(y=ratioenrichment,x=reorder(pathway,as.numeric(ratioenrichment))))+geom_bar(stat="identity",aes(fill=source))+coord_flip()+
  xlab("")+ylab("Ratio of enrichment")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))+ggtitle("GFP-Out/GFP-In")+scale_fill_manual(values=c("royalblue4","lightseagreen"))

