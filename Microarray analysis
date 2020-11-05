library(ggplot2)
library(reshape2)
library(limma)

arguements = commandArgs(T)
setwd(arguements[1])
rd<-read.csv(arguements[2], header = T)
genenames<-rd[,1]
rd<-rd[,-1]

png(filename = paste(arguements[2],"without_norm.boxplot.png"), width = 800, height = 800)
boxplot(rd)
dev.off()

rd<-backgroundCorrect(rd, method="normexp", offset=0)
rd<-normalizeBetweenArrays(y, method="quantile")

png(filename = paste(arguements[2],"after_norm.boxplot.png"), width = 800, height = 800)
boxplot(y)
dev.off()


a<-apply(rd,2,as.numeric)
rownames(a)<-genenames
aaverage<-avereps(as.matrix(a))
Tumor <- rep("T",times=arguements[3])
Control <- rep("N",times=arguements[4])
treatment <- c(Tumor,Control)
treatment1 <- factor(treatment)
design1 <- model.matrix(~0 + treatment1)
colnames(design1) <- c("N","T")
cntrast.matrix1 <- makeContrasts(N-T,levels = design1)
fit1 <- lmFit(aaverage,design1)
fit1 <- contrasts.fit(fit1, cntrast.matrix1)
fit12 <- eBayes(fit1)
dif1 <- topTable(fit12,coef = 1,n=nrow(fit12))
write.csv(dif1, paste(arguements[2],'DE.csv'))

dif1$threshold = factor(ifelse(dif1$adj.P.Val < 0.05 & abs(dif1$logFC) >= 1, ifelse(dif1$logFC>= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(dif1,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
  geom_point(alpha=0.4, size=2.5)+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  theme_bw()+
  theme(legend.title = element_blank())+
  ylab('-log10 (adj.Pval)')+xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)

library('pheatmap')
degs<-read.csv('DE gene exp.csv', header = T)
degs<-aggregate(.~Gene.symbol, degs, mean)
deg_names<-degs[,1]
degs<-degs[,-1]
rownames(degs)<-deg_names
annotation_col = data.frame(Type = factor(rep(c("Tumor", "Normal"), c(52, 6))), Samples = colnames(degs))
rownames(annotation_col) = paste(colnames(degs))
png(filename = paste('DEGs',"heatmap.png"), width = 1500, height = 2000)
pheatmap(degs, scale = "row", clustering_method = "average", clustering_distance_rows = "correlation", color = colorRampPalette(c("navy", "white", "firebrick3"))(50), show_rownames=F, show_colnames = F, annotation_col = annotation_col)
dev.off()

