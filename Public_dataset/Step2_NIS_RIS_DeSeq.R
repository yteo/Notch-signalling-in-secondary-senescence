
#Differentially expressed genes in NIS-RIS

library(DESeq2)

directory = ">"

# NIS vs RIS
sampleFiles<-grep('.*_sorted_htseq.txt ',list.files(directory),value=TRUE)
#countFile = paste0(sample_name, "_sorted_htseq.txt")
condition<-factor(sub("\\_.*","", sampleFiles))
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = condition
                          
)


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = directory,
                                          design = ~ condition)

DESeq2Table <- DESeq2Table[ rowSums(counts(DESeq2Table)) > 1, ]
DESeq2Table$condition <- factor(DESeq2Table$condition, levels=c("RIS","NIS"))


dds<-DESeq2Table
dds <- estimateSizeFactors(dds)
ddscount<-counts(dds, normalized=TRUE)

DESeq2Table <- DESeq(DESeq2Table)
res <- results(DESeq2Table)
res<-res[order(res$padj),]

res2<-data.frame(gene_id=rownames(res),log2FC=res$log2FoldChange,padj=res$padj)
gene_name<-read.table("AllRNA_ID_name",header=T)
res2<-merge(gene_name,res2)
write.csv(res2,file="NIS_RIS.csv",quote=F)
