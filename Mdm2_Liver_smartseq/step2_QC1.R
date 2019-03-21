library(ggplot2)

df<-read.table("mappingstat2.txt",header=T,row.names=1)
#df<-df[grep("6_1",rownames(df)),]
allcond<-sub("6_1b.*","6_1b",sub("6_1c.*","6_1c",sub("6_2b.*","6_2b",sub("6_2c.*","6_2c",rownames(df)))))
df$con<-allcond
df$cri<-rep("failedERCCratio",length(allcond))
# passed ERCC ratio if ERCC mapped / total mapped <0.5
df[which((df$chrmapped>50000 )),]$cri <- "passed"


hist(log10(df$chrmapped),breaks=200,ylab="Number of Mdm2 deleted cells",xlab="log mapped reads",main="Number of mapped reads")
abline(v=log10(50000),lty=2,col="red")
text(5.5,6,"75 cells",col="red")

passed<-df[which(df$cri=="passed"),]
write.csv(rownames(passed) ,file="50000_mappedreads.txt",quote=F,row.names = F)


