
library(ggplot2)

df<-read.table("mappedreadsstat.txt",header=F,row.names=1)
allcond<-sub("Sen1","Sen",sub("Sen3","Sen",sub("Sen2","Sen",sub("Growing2","Growing",gsub("\\_.*","",rownames(df))))))
df$con<-allcond
df$cri<-rep("failedERCCratio",length(allcond))
# passed ERCC ratio if ERCC mapped / total mapped <0.5
df[which((df$V3 < (0.5*df$V1) & df$V2>200000 & df$V3>0)),]$cri <- "passedERCCratio"



#########################################
#Plot ERCC/total mapped reads ratio versus chromosmoe mapped reads with condition


ggplot()+geom_point(aes(y=(df$V3/df$V1),x=log10(df$V2),color=factor(allcond)))+theme_bw()+xlab("log 10 Chromosome mapped reads")+
  ylab("ERCC ratio")+geom_vline(xintercept=log10(200000),color="red")+
  geom_hline(yintercept=0.5,color="red")


########################################

#Plot total mapped reads versus chromosome mapped reads. colored by 
#1. chromosome mapped reads > 200,000
#2. ERCC/total mapped reads ratio <0.5

ggplot()+geom_point(aes(y=(log10(df$V1)),x=log10(df$V2),color=factor(df$cri)))+theme_bw()+xlab("log 10 Chromosome mapped reads")+
  ylab("log 10 total mapped")+geom_vline(xintercept=log10(200000),color="red")
#########################################
#Plot ERCC mapped read versus chromosome mapped reads. colored by

#1. chromosome mapped reads > 200,000
#2. ERCC/total mapped reads ratio <0.5
ggplot()+geom_point(aes(y=log10(df$V3),x=log10(df$V2),color=factor(df$cri)))+theme_bw()+xlab("log10 Chromosome mapped reads")+
  ylab("log10 ERCC mapped reads")+geom_vline(xintercept=log10(200000),color="red")

#########################################
#Plot the total hg19 mapped reads for each condition

logy<-log10(df$V2)
ggplot(df,aes(x=reorder(rownames(df),-logy),y=logy))+geom_bar(aes(fill=con),stat="identity")+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
  ylab("Log10 (Total chr mapped)")+xlab("Conditions")+
  ggtitle("Total chr mapped in scRNA samples in log10 scale")+
  geom_hline(yintercept=log10(200000),linetype="dashed",color="black",size=1)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),panel.background=element_blank(),panel.border=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  scale_fill_hue(c=80,l=80)+scale_fill_brewer(palette="Spectral")

###########################################
#Cells passed requirements:
#a) chromosome mapped reads > 200,000
#b) ERCC/total mapped reads ratio <0.5

passed<-df[which(df$cri=="passedERCCratio"),]

write.csv(passed,file="passedmappedreads.csv",rownames=F,quote=F)
