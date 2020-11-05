library(edgeR)
library(stringi)
library(data.table)

rm(list=ls())
miRNA_exp<-read.csv('exp_matrix_grouped.csv', header = T)

##放弃merge匹配了，直接excelvlookup匹配好生成对应的groupinfo.csv文件，继续
#grp<-as.factor(merge(exp_mat_file_list, cases_annot[, c("file_id", "cases.0.case_id", "cases.0.samples.0.sample_type")], by="file_id")[,cases.0.samples.0.sample_type])

groupinfo<-read.csv('group_info.csv', header = T)
grp<-as.factor(groupinfo[,"sample_type"])
table(grp)

miRNA_names<-miRNA_exp[,1]
miRNA_exp<-miRNA_exp[,-1]
rownames(miRNA_exp)<-miRNA_names
miRNA_exp_2<-copy(miRNA_exp)
grp_chr<-as.character(grp)
grp_chr[which(grp_chr=="Primary Tumor")]<-"T"
grp_chr[which(grp_chr=="Solid Tissue Normal")]<-"N"
grp<-as.factor(grp_chr)
colnames(miRNA_exp_2)<-paste(grp,1:length(grp),sep="")
pairs(log2(1+miRNA_exp_2[,1:7]), pch=".",lower.panel=NULL) #查看前7个样本表达情况，样本间的表达量基本都是呈线性的。

d <- DGEList(counts=miRNA_exp_2, group=grp)
d <- calcNormFactors(d,method = "TMM")
head(d$samples)

cols <- as.numeric(d$samples$group)
plotMDS(d,col=cols) #离群样本最好去除，标准情况是N和N聚在一块，T和T聚在一块，要把两个聚类中对方的样本去除。

mm <- model.matrix(~-1+grp)
d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)
sqrt(d$common.dispersion) #该值最好在0.2 - 0.4之间。
plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE) #可能会报错
plotSmear(d, pair=c("N","T"), ylim=c(-5,5))
f <- glmFit(d,mm)
con <- makeContrasts("T-N"=grpT-grpN,levels=colnames(mm))
lrt <- glmLRT(f,contrast=con)
miRNA_DE_result<-topTags(lrt,Inf)$table
miRNA_DE_result<-cbind(miRNAs=rownames(miRNA_DE_result),miRNA_DE_result)
write.csv(miRNA_DE_result, 'edgeR_miRNA_DE.csv')
