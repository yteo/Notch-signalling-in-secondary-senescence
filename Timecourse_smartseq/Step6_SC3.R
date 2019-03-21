library(sc3)
# split cells is the file that contains sample name for senbottom and sentop
splitcells<-read.csv("split.csv",header=T)
# Sen_17cells_individualbranch is senescent botttom branch from monocle (secondary senescence)
# Sen_22cells_day24branch is senescenct top branch from monocle(OIS)

dat<-read.csv("filtered_raw_IMR90.txt")
gene_name<-read.table("AllRNA_ID_name",header=T)
dat<-merge(gene_name,dat,by.y="row.names",by.x="gene_id")
rownames(dat)<-dat[,2]
rownames(dat)<-make.names(dat$name,unique=T)
dat<-dat[,-1]
dat<-dat[,-1]
dat<-data.frame(dat)
dat<-dat[,grep("Sen",colnames(dat))]


SenBottom<-dat[, colnames(dat) %in% splitcells$Sen_17cells_individualbranch]
colnames(SenBottom) <- paste("SenBottom", colnames(SenBottom), sep = "_")

Senday24<-dat[, colnames(dat) %in% splitcells$Sen_22cells_day24branch]
colnames(Senday24) <- paste("SenTop", colnames(Senday24), sep = "_")

others<-dat[,!(colnames(dat) %in% splitcells$Sen_17cells_individualbranch)]
others<-others[, !(colnames(others) %in% splitcells$Sen_22cells_day24branch)]
colnames(others) <- paste("NA", colnames(others), sep = "_")

dat_split<-cbind(SenBottom,Senday24,others)

dat_split_cell_info$top_bottom<-factor(gsub("(SenTop|SenBottom|NA).*", "\\1", colnames(dat_split)))



dat_split_sceset@sc3 <- list()
dat_split_cell_exprs <- as.matrix(sapply(dat_split, as.numeric)) 
rownames(dat_split_cell_exprs)<-rownames(dat_split)

pd <- new("AnnotatedDataFrame", data = dat_split_cell_info)
dat_split_sceset <- newSCESet(countData = dat_split_cell_exprs, phenoData = pd)

# vary the cluster numbers to see which one is the best. 
dat_split_sceset<-sc3(dat_split_sceset, ks =2)