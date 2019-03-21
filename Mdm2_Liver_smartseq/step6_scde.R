# scde of induced versus control (WT) 
# continuation from step5_QC2_boxplot_heatmap.R
dat_MDM<-dat[,colnames(dat) %in% gsub(".*6_","6_",colnames(exprs_mat2))]
dat_WT<-dat[,colnames(dat) %in% gsub(".*6_","6_",colnames(exprs_mat3))]
colnames(dat_MDM)<-paste("MDM_",colnames(dat_MDM),sep="")
colnames(dat_WT)<-paste("WT_",colnames(dat_WT),sep="")
mdm_wt<-cbind(dat_MDM,dat_WT)

cd_sc3<-mdm_wt
# omit genes that are not expressed in sufficient number of cells (more than 50)
cd_sc3 <- cd_sc3[rowSums(cd_sc3>0)>2,]

colnames(cd_sc3)<-gsub("X","",colnames(cd_sc3))

con_filter <- factor(gsub("(MDM|WT).*", "\\1", colnames(cd_sc3)), levels = c("MDM", "WT"))
names(con_filter) <- colnames(cd_sc3)  
table(con_filter)


# calculate models


o.ifm <- scde.error.models(counts = cd_sc3, groups = con_filter, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)


# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd_sc3, length.out = 400, show.plot = FALSE)

# define two groups of cell
groups <- factor(gsub("_.*","", rownames(o.ifm)), levels  =  c("MDM","WT"))
#groups <- factor(gsub("(1|3_4).*", "\\1", rownames(o.ifm)), levels  =  c(paste(cluster1),paste(cluster2))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd_sc3, o.prior,groups=groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
significant.genes <- which(p.values.adj<1.1)
length(significant.genes)


ord <- order(p.values.adj[significant.genes]) # order by p-value
de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
colnames(de) <- c("Lower bound","log2FC","Upper bound","p-value")
ensembl<-read.table("modified_ENSM_genename_mm10.txt")
dat2<-merge(de,ensembl,by.x="row.names",by.y="V1")

write.table(dat2,file="SCDE_MDM_WT.txt",quote=F,sep="\t")
write.table(ediff,file="SCDEEDIFF_MDM_WT.txt",quote=F,sep="\t")
write.table(dat2$V2,file="mdm_wt_background.txt",quote=F,row.names = F)



########## 
#webgestalt
# genes with logFC>0 , padj < 0.05 between MDM exp / WT
pathway<-data.frame(pathway=c("p53 signaling pathway","Axon guidance","Cytokine-cytokine receptor interaction","miRNA regulation of DNA Damage Response","p53 signaling"), source=c("Kegg","Kegg","Kegg","WikiPathway","WikiPathway"),ratioenrichment=c("5.31","3.49","3.78","5.07","4.96"))


library(ggplot2)
ggplot(data=pathway,aes(y=ratioenrichment,x=reorder(pathway,as.numeric(ratioenrichment))))+geom_bar(stat="identity",aes(fill=source,color=source))+coord_flip()+
  xlab("")+ylab("Ratio of enrichment")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))



# 
expr.count_passed<-expr.count[,colnames(expr.count) %in% passed]
write.table(expr.count_passed,file="Hepatocyte_htseq.txt",sep="\t",quote=F)







#############################################################################
#SCDE of secondary and primary senescence in induced experiment

cd_sc3<-exprs_mat2

cd_sc3<-cd_sc3[rowSums(cd_sc3)>0, ] 
cd_sc3<-data.frame(cd_sc3)


cd_sc3 <- cd_sc3[rowSums(cd_sc3>0)>2,]

colnames(cd_sc3)<-gsub("X","",colnames(cd_sc3))

con_filter <- factor(gsub("(MdmPos|MdmNeg).*", "\\1", colnames(cd_sc3)), levels = c("MdmPos", "MdmNeg"))
names(con_filter) <- colnames(cd_sc3)  
table(con_filter)


# calculate models


  o.ifm <- scde.error.models(counts = cd_sc3, groups = con_filter, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

save(o.ifm,file="./SCDE_byexon/o.ifm.Rdata")



# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd_sc3, length.out = 400, show.plot = FALSE)

# define two groups of cell
groups <- factor(gsub("_.*","", rownames(o.ifm)), levels  =  c("MdmPos","MdmNeg"))
#groups <- factor(gsub("(1|3_4).*", "\\1", rownames(o.ifm)), levels  =  c(paste(cluster1),paste(cluster2))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd_sc3, o.prior,groups=groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)


# convert z score to p value
scde.test.gene.expression.difference("Maml1", models = o.ifm, counts = cd_sc3, prior = o.prior)
scde.test.gene.expression.difference("Rfng", models = o.ifm, counts = cd_sc3, prior = o.prior)
scde.test.gene.expression.difference("Smad3", models = o.ifm, counts = cd_sc3, prior = o.prior)
#http://hms-dbmi.github.io/scw/differential-expression.html
p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
significant.genes <- which(p.values.adj<1.1)
length(significant.genes)


ord <- order(p.values.adj[significant.genes]) # order by p-value
de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
colnames(de) <- c("Lower bound","log2FC","Upper bound","p-value")
de<-de[order(de$log2FC),]

## mouse and human
hm<-read.table("human_mouse_genename.txt",header=T,sep="\t")
de_hm<-merge(hm,de,by.x="Gene.name",by.y="row.names")

gsea<-data.frame(gene=de_hm$Human.gene.name,log2FC=de_hm$`log2FC`)

#######

p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR

ediff<-cbind(ediff[,1:3],ediff[,5],ediff[,6],p.values.adj)
colnames(ediff) <- c("Lowerbound","log2FC","Upperbound","Zscore","Corrected_Zscore","adj.p")

#### only mouse gene name



gsea<-data.frame(gene=rownames(ediff),log2FC=ediff$`log2FC`)
write.table(gsea,file="./SCDE_byexon/MdmPosneg_exon56_log2FC_final_mouse.rnk",quote=F,sep="\t",row.names=F)




write.table(ediff,file="./SCDE_byexon/SCDE_mdm2sepbyExonPOS_NEG_final_mouse.txt",quote=F,sep="\t",row.names=T)




##################
# webgestalt
# ratio of enrichment

path=data.frame( wikipath=c("Notch Signaling Pathway","ESC Pluripotency Pathways","Delta-Notch Signaling Pathway","Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling","miRNA regulation of DNA Damage Response","Nuclear Receptors","IL-2 Signaling Pathway","TGF Beta Signaling Pathway","ErbB signaling pathway","PluriNetWork"),ratioenrichment=
c(6.39,4.34,4.19,6.77,3.92,5,3.03,3.71,3.6,1.77))


#min:20 genes
#max: 400 genes
path=data.frame( source=rep("Wikipathway",3),pathway=c("Notch Signaling Pathway","Delta-Notch Signaling Pathway","TGF Beta Signaling Pathway"),ratioenrichment=
                   c(7.07,4.63,4.11))

kegg=data.frame(source="Kegg",pathway=c("Notch signaling pathway"),ratioenrichment="5.21")

path<-data.frame(rbind(path,kegg))
library(ggplot2)
ggplot(data=path,aes(y=ratioenrichment,x=reorder(pathway,as.numeric(ratioenrichment))))+geom_bar(stat="identity",aes(fill=source,color=source))+coord_flip()+
  xlab("")+ylab("Ratio of enrichment")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(colour = "black"),
axis.title.y = element_text(colour = "black"))