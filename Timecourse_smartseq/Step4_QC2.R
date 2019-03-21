### IMR 90 Filtering based on gene counts ####

library("ggplot2")
library(matrixStats)
library(tidyr)
dat<-read.table("Compiled_downsampled-chr.txt")
dat<-spread(dat,V1,V3)
rownames(dat)<-dat$V2
dat<-dat[,-1]



rownames(dat)<-sub("*\\..*", "", rownames(dat))


expr.count <- dat[grep("^ERCC",row.names(dat),invert=TRUE),]

#Sum the count in the ERCC controls and calculate the ratio with the total reads in mRNAs

total.sum.counts <- apply(dat,2,sum)



#Calculate the number of genes with at least 1 read in each sample (cell)


genes.with.reads <- sapply(1:ncol(expr.count),function(x){sum (!expr.count[,x] == 0)})
names(genes.with.reads) <- colnames(expr.count)


#Plot the total read count in mRNAs to find how many cells have low read count

options(scipen=999)

par(mfrow = c(2,2))
hist(total.sum.counts,breaks = 100, ylim=c(0,200))
hist(total.sum.counts,breaks = 200,xlim = c(0,50000000),ylim=c(0,200))
hist(total.sum.counts,breaks = 1000,xlim = c(0,10000000),ylim=c(0,200))
hist(total.sum.counts,breaks = 2000,xlim = c(0,50000000),ylim=c(0,200))
par(mfrow = c(1,1))

hist(total.sum.counts[total.sum.counts > 100000],breaks = 100,ylim=c(0,200))

options(scipen=-1)


#Plot for each cell the number of genes showing at least 1 read

plot(total.sum.counts,genes.with.reads,pch=20,xlab="Total Read Count",ylab="Genes with at least 1 read")



#Set up the quality filters


filter.ERCC <- as.logical(ratio.ERCC >0 & ratio.ERCC <0.5 )
filter.genecounts <- as.logical(genes.with.reads > 500)
filter.total.reads <- as.logical(total.sum.counts > 100000)


#Plot num genes with at least 1 read vs total gene count

############################
con <- factor(gsub("(Sen|day2|day4|Growing).*", "\\1", colnames(dat)), levels = c("Sen","day2","day4","Growing"))

ggplot()+geom_point(aes(x=total.sum.counts,y=genes.with.reads,color=con))+xlab("Total gene count")+ylab("genes with at least one read")+geom_vline(xintercept=80000,linetype="dashed")+geom_hline(yintercept=500,linetype="dashed")+theme_bw()+scale_color_manual(values=c("red","palegreen4","plum4","skyblue"))+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    panel.background = element_blank()) 

#Cells that passed filter

name<-as.data.frame(total.sum.counts)
name$con<-row.names(name)
passed<-name[( genes.with.reads > 500 & total.sum.counts > 80000),]$con
write.csv(passed,file="passedgenecount.csv")

#compile gene counts for samples that passed filter, name as "Compiled_downsampled_chr.txt"

passed_dat<-dat[,colnames(dat) %in% passed]
write.table(passed_dat,file="Compiled_downsampled_chr.txt",quote=F,sep="\t")