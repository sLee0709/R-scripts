library(ROCR)
library(ggpubr) 
library(patchwork)
library(caret)
library(glmnet)
library(ggfortify)
library(pheatmap)
library(factoextra)
library(FactoMineR)

rm(list=ls())

#读入数据
exprSet<-read.csv('candidate_biomarkers.csv', header = T)
meta<-read.csv('survival_info.csv', header = T)

rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-1]

set.seed(2)
sam<- createDataPartition(meta$status_code, p = .7,list = FALSE)          #这里的p是指训练集百分比，list=F是不转换为list

#训练集和测试集划分
train <- exprSet[,sam]
test <- exprSet[,-sam]
train_meta <- meta[sam,]
test_meta <- meta[-sam,]

x = t(log2(train+1))
y = train_meta$status_code
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 50,alpha = 1)       #这里的nlambda是指让算法自动挑选n个不同的λ值，拟合出n个系数不同的模型。
plot(cv_fit)

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)   #选取lambda.min建模
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)   #选取lambda.1se建模

lasso.prob <- predict(cv_fit, newx=t(log2(test+1)), s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(test_meta$status_code ,lasso.prob)
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)

#查看模型对测试集的区分度
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()

#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")

tpr_min = performance(pred_min,"tpr")@y.values[[1]]
tpr_1se = performance(pred_1se,"tpr")@y.values[[1]]
dat = data.frame(tpr_min = perf_min@y.values[[1]],
                 fpr_min = perf_min@x.values[[1]],
                 tpr_1se = perf_1se@y.values[[1]],
                 fpr_1se = perf_1se@x.values[[1]])
ggplot() + 
  #geom_line(data = dat,aes(x = fpr_min, y = tpr_min),color = "blue") + 
  geom_line(data = dat,aes(x = fpr_1se, y = tpr_1se),color = "red")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  #annotate("text",x = .75, y = .25, label = paste("AUC of min = ",round(auc_min,2)),color = "blue")+
  annotate("text",x = .75, y = .15,label = paste("AUC of 1se = ",round(auc_1se,2)),color = "red")+
  scale_x_continuous(name  = "fpr")+
  scale_y_continuous(name = "tpr")

fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]

#读取全部的miRNA和样本matrix
expr<-read.csv('liqun_removed_matrix.csv', header = T)
rownames(expr)<-expr[,1]
expr<-expr[,-1]
tumor<-rep('tumor',529)
normal<-rep('normal', 71)
group_list<-c(tumor, normal)

choose_matrix=expr[choose_gene,]
n=t(scale(t(log2(choose_matrix+1))))
n[n>2]=2
n[n< -2]= -2
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)

#特征基因热图
pheatmap(n,show_colnames = F,annotation_col = annotation_col,
         filename = 'lasso_genes.heatmap.png', clustering_method = 'average')

#PCA
df=as.data.frame(t(choose_matrix))
df$group=group_list
df.pca=PCA(df[,1:(ncol(df)-1)], graph = FALSE)
png('lasso_genes.pca.png',res=120)
fviz_pca_ind(df.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df$group, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()
