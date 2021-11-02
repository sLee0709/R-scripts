"

通过TCGA表达矩阵获取DE基因
开始前准备：1.表达矩阵；2.样本信息（两列表格，第一列为样本名称，第二列为分类标签，如Tumor，Normal等等。）

"

library(tidyverse)
library(DESeq2)
mycounts<-read.csv('mycounts.csv', header = T)
row.names(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
conditions<-read.csv('mycounts_condition.csv', header = F)
conditions$V2<-factor(conditions$V2, c('Tumor','Normal'))
dds<-DESeqDataSetFromMatrix(countData = mycounts, colData = conditions, design = ~V2)
dds<-DESeq(dds)

"
运行log ↓↓↓

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 2183 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing

"

sizefactor<-sizeFactors(dds)
resultsNames(dds)
res<-results(dds,name="V2_Normal_vs_Tumor")
head(res)
res<-res[order(res$padj),]
head(res)
resDF<-as.data.frame(res)
resDF$Gene_id=row.names(resDF)
head(resDF)
resDF<-resDF[,c(7,1:6)]
head(resDF)
write.csv(resDF, 'protein_coding_DE.csv')
