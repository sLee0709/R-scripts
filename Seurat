setwd('D:/WenLabTask')

library(SeuratDisk)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)
library(clustree)
library(ggpubr)
library(tidyverse)
library(ggstatsplot)

#Convert('TS_Vasculature.h5ad','h5seurat',overwrite = F)
seurat_obj<-LoadH5Seurat('TS_Vasculature.h5seurat', meta.data = T, misc = F)
#seurat_obj<-subset(seurat_obj_all, subset = (donor == 'TSP1' | donor == 'TSP3' | donor == 'TSP4' | donor == 'TSP5' | donor == 'TSP6' | donor == 'TSP7' | donor == 'TSP8' | donor == 'TSP9' | donor == 'TSP10' | donor == 'TSP12' | donor == 'TSP15')) #样本量过大，从数据集中挑几个sample出来分析。
#seurat_obj@meta.data<-separate(data = seurat_obj@meta.data, col = development_stage, into = c('age','age0'), sep = '-y')
#seurat_obj@meta.data<-seurat_obj@meta.data[,-29]
#seurat_obj@meta.data$age_class<-ifelse(seurat_obj@meta.data$age > 40, ifelse(seurat_obj@meta.data$age > 60,'Old','Mid-age'),'Young')
#or
#seurat_obj@meta.data$age_class<-ifelse(seurat_obj@meta.data$donor == 'TSP10' | seurat_obj@meta.data$donor == 'TSP9' | seurat_obj@meta.data$donor == 'TSP4' | seurat_obj@meta.data$donor == 'TSP5','Young',ifelse(seurat_obj@meta.data$donor == 'TSP2' | seurat_obj@meta.data$donor == 'TSP6' | seurat_obj@meta.data$donor == 'TSP7' | seurat_obj@meta.data$donor == 'TSP12','Old','Mid-age'))

#查看是否存在批次效应
#DimPlot(seurat_obj, reduction = 'umap', group.by = 'donor')

#View(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)[3]<-"nCount_RNA"
colnames(seurat_obj@meta.data)[4]<-"nFeature_RNA"
seurat_obj[['percent.mt']]<-PercentageFeatureSet(seurat_obj, pattern = "^MT-") #对于去除batch effect后的integrated数据报错，原因和解决方案：https://github.com/satijalab/seurat/issues/1665
seurat_obj<-subset(seurat_obj, subset = nCount_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 30 & nFeature_RNA > 600 )
VlnPlot(seurat_obj,features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seurat_obj<-NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
seurat_obj<-FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)

#去除批次效应（if applicable）
#See batch_effect_removement.R


top10<-head(VariableFeatures(seurat_obj), 10)
plot1<-VariableFeaturePlot(seurat_obj)
plot2<-LabelPoints(plot = plot1, points = top10, repel = F)
plot2

all.genes<-rownames(seurat_obj)
seurat_obj<-ScaleData(seurat_obj, features = all.genes) #如果报错空间不足xxxGB，执行memory.limit(size=56000)
#seurat_obj<-ScaleData(seurat_obj, vars.to.regress = 'percent.mt')
seurat_obj<-RunPCA(object = seurat_obj, pc.genes=VariableFeatures(seurat_obj))
seurat_obj<-FindNeighbors(seurat_obj,dims = 1:10)

#执行FindClusters之前用clustree看一下resolution的最优选择（即细胞集分多成少个cluster最合适）
#clustfind<-FindClusters(seurat_obj, resolution = seq(0.5, 1.2, by=0.1))
#clustree(clustfind)
"
label_position<-function(labels){
  if(length(unique(labels)) == 1){
    position<-as.character(unique(labels))
  }else {
    position<-'mixed'
  }
  return(position)
}

"
#clustree(clustfind@meta.data, prefix = 'RNA_snn_res.', node_label = 'cell_ontology_class', node_label_aggr = "label_position") #这里cell_ontology_class是先验的细胞类型，如果没有需要先给细胞添加注释。

seurat_obj<-FindClusters(seurat_obj, resolution = 0.5)
seurat_obj<-RunUMAP(seurat_obj, dim=1:10)
DimPlot(seurat_obj, reduction = 'umap', label = T)

#细胞注释
#这个数据集已经分好细胞类型，所以group.by直接填写'cell_ontology_class'即可，否则需要执行下面代码进行细胞注释
obj_singleR<-GetAssayData(seurat_obj, slot='data')
hpca.se<-HumanPrimaryCellAtlasData()
cellType_map<-SingleR(test = obj_singleR, ref = hpca.se, labels = hpca.se$label.main)
seurat_obj@meta.data$labels<-cellType_map$labels
DimPlot(seurat_obj, group.by = c("seurat_clusters", "labels"),reduction = "umap", label = T)

#查看endothelin相关基因表达情况。
endogene<-c('EDN1', 'EDNRA', 'EDNRB', 'ACE2', 'ACTA2', 'CDKN2A', 'CDKN1A', 'TP53', 'MDM2')
"
png('MDM2.png', units = 'in', width = 15, height = 10, res = 300)
VlnPlot(seurat_obj, features = 'MDM2', group.by = 'cell_ontology_class', split.by = 'age_class',cols = c('red','blue','green'))
dev.off()

png('FeaturePlot.png', units = 'in', width = 15, height = 15, res = 300)
FeaturePlot(seurat_obj, features = endogene, reduction = 'umap', pt.size = 2)
dev.off()

"
VlnPlot(seurat_obj, features = endogene, group.by = 'cell_ontology_class')
FeaturePlot(seurat_obj, features = endogene, reduction = 'umap', pt.size = 1)
DotPlot(seurat_obj, features = endogene, group.by = "cell_ontology_class") + coord_flip() + RotatedAxis()

#比较young和old的cell分布
seurat_obj@meta.data<-unite(seurat_obj@meta.data, "age_cell", age_class, cell_ontology_class, sep = '-', remove = F)
Young<-WhichCells(seurat_obj, cells = grep('Young-vein endothelial cell', seurat_obj$age_cell, value = F))
Midage<-WhichCells(seurat_obj, cells = grep('Mid-age-vein endothelial cell', seurat_obj$age_cell, value = F))
Old<-WhichCells(seurat_obj, cells = grep('Old-vein endothelial cell', seurat_obj$age_cell, value = F))

"
png('FeaturePlot.png', units = 'in', width = 5, height = 5, res = 300)
DimPlot(seurat_obj,reduction = 'umap', label = F, cells.highlight = list(Young, Midage, Old)) + scale_color_manual(labels = c('Other', 'Old', 'Mid-age', 'Young'), values = c('grey', 'red', 'blue', 'green')) + labs(title = 'endothelial cell of artery') + theme(plot.title = element_text(hjust = 0.5))
dev.off()
"

DimPlot(seurat_obj,reduction = "umap", label = F, cells.highlight = list(Young, Midage, Old)) + scale_color_manual(labels = c("Other", "Old", "Mid-age", "Young"), values = c("grey", "red", 'blue', 'green')) + labs(title = 'endothelial cell of artery') + theme(plot.title = element_text(hjust = 0.5))

#计算某一类细胞Young和Old的差异表达基因
Idents(cap_obj2)<-'age_class'
cap_exp<-FindMarkers(cap_obj2, ident.1 = 'Young', ident.2 = 'Old', verbose = F)

#可视化某特定基因在某细胞中的表达情况（boxplot）
vp_case1 <- function(dataset, gene_signature, test_sign, y_max, c1='red', c2='blue', c3='green'){
  plot_case1 <- function(signature){
    VlnPlot(dataset, features = signature,
            pt.size = 0.1, 
            group.by = "age_factor", 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c(c1, c2, c3)
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif") 
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
}

gene_sig<-'EDNRB'
comparison<-list(c('Young','Old'),c('Young','MidAge'),c('MidAge','Old'))
cap_obj@meta.data$age_factor<-factor(x = cap_obj@meta.data$age_class, levels = c('Young', 'MidAge', 'Old')) # To adjust the order of violin plot.
vp_case1(dataset = cap_obj, gene_signature = gene_sig, test_sign = comparison, y_max = 6, c1 = '#8FBC8F',c2 = '#F0E68C', c3 = '#F08080')

#gene-gene相关性
gene_sig<-c('EDN1', 'EDNRA', 'EDNRB','SIRT1','SIRT6','SIRT7')
exp_mt<-seurat_obj@assays$RNA@data
exp_mt<-as.data.frame(exp_mt)
exp_mt<-t(exp_mt)
exp_mt<-as.data.frame(exp_mt)
y<-as.numeric(exp_mt[,'EDNRB'])
colname<-colnames(exp_mt)
cor_data_df<-data.frame(colname)
for (i in 1:length(colname)) {
  test<-cor.test(as.numeric(exp_mt[,i]),y,type = 'spearman')
  cor_data_df[i,2]<-test$estimate
  cor_data_df[i,3]<-test$p.value
}
names(cor_data_df)<-c('Gene','correlation','pvalue')
write.csv(cor_data_df,'correlations.EDNRB-allgene.csv')

##调取top10正负相关的基因
top10_negative<-cor_data_df[order(cor_data_df[,2]),][1:10,]
top10_positive<-cor_data_df[order(-cor_data_df[,2]),][1:11,]
top10_positive<-top10_positive[-1,]
top10<-rbind(top10_positive,top10_negative)
##作图
##（顺便先看一下SIRT1，6，7的表达情况）
SIRT<-exp_mt[,c('SIRT1','SIRT6', 'SIRT7','EDNRA')]
png('EDNRB-SIRT.exp.png', units = 'in', width = 5, height = 4, res = 300)
boxplot(SIRT)
dev.off()

"
#将EDN1,EDNRA,EDNRB和SIRT1/6/7中含有0值的行删除后做相关性分析
exp_mt2<-exp_mt[,gene_sig]
exp_mt3<-exp_mt2[exp_mt2$EDNRB>0,]
exp_mt3<-exp_mt3[exp_mt3$SIRT7>0,]
png('EDNRB-SIRT7.corr.zero-removed.png', units = 'in', width = 4, height = 4, res = 300)
ggscatter(exp_mt3, x='EDNRB', y='SIRT7',
          color = 'black', size = 2,
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgrey'),
          conf.int = T,
          cor.coef = T,
          cor.coeff.args = list(method = 'spearman', label.x = 0.5, label.sep = '\n'))
dev.off()

"

png('EDNRB-SIRT7.corr.png', units = 'in', width = 4, height = 4, res = 300)
ggscatter(exp_mt, x='EDNRB', y='SIRT7',
          color = 'black', size = 2,
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgrey'),
          conf.int = T,
          cor.coef = T,
          cor.coeff.args = list(method = 'spearman', label.x = 1, label.sep = '\n'))
dev.off()

top10$Group<-ifelse(top10$correlation >0, 'Positive','Negative')
png('t10.corr.genes.png', units = 'in', width = 5, height = 5, res = 300)
ggplot(top10, aes(x = correlation,y = reorder(Gene, correlation))) + 
  geom_bar(stat = 'identity', aes(fill=Group)) + 
  theme_bw() + 
  ylab('Gene Symbols') + labs(fill='Correlation') + 
  scale_fill_manual(values = c('#5F9EA0','#FFA07A'))
dev.off()

#调取cluster的feature基因
cap_c1<-subset(cap_obj, subset = seurat_clusters == '1')
cap_c1<-FindVariableFeatures(cap_c1, selection.method = 'vst', nfeatures = 100)
t100<-head(VariableFeatures(cap_c1), 100)
write.csv(t100, 'cap.c1.features.csv')
rm(cap_c2, t100)

findfeatures<-function(clst.idents, ft.nums, top.nums){
  clst<-subset(cap_obj, subset = seurat_clusters == clst.idents)
  clst<-FindVariableFeatures(clst, selection.method = 'vst', nfeatures = ft.nums)
  tops<-head(VariableFeatures(clst), top.nums)
  write.csv(tops, 'cap.clst.features.csv')
}

findfeatures(clst.idents = 6, ft.nums = 100, top.nums = 100)



