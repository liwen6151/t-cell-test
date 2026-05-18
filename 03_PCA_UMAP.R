rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("03_PCA_UMAP")) {dir.create("03_PCA_UMAP")}
setwd("03_PCA_UMAP/")


#PCA UMAP降维聚类-------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(scales)
#library(harmony)
library(ggplot2)
#library(GSEABase)
library(tidyr)

load("../02_VariableFeatures/scRNA_VariableFeatures.Rdata")

#在上一步已经做了scaledata
scRNA.nor.sca <- scRNA
scRNA.norm.pca <- RunPCA(scRNA.nor.sca, features = VariableFeatures(object = scRNA.nor.sca), npcs = 50)

print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
##可视化对每个成分影响比较大的基因集,自己了解就行，1：2~50都可以
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")


library(tidyr)

pdf("02.PCA.pdf",w=6,h=5,family = 'Times')
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'orig.ident')
dev.off()
png("02.PCA.png",w=6,h=5,units = 'in',res = 600,family='Times')
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'orig.ident')
dev.off()


## 首先找到最佳聚类数

## JackStraw函数的功能/作用：
## 随机排列数据的子集，并计算这些“随机”基因的预测PCA分数。
## 然后将“随机”基因的PCA分数与观察到的PCA分数进行比较，以确定统计显著性。
## 最终结果是每个基因与每个主成分相关的p值
## 简而言之: JackStraw()函数可以计算出每个主成分中各基因的P值，用于判断哪些主成分更具有统计学意义
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) # num.replicate : 要执行的复制采样数；dims : 计算重要性的PC数
## ScoreJackStraw()函数：用于量化主成分的显著性强度，富含低P值基因较多的主成分更有统计学意义
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
system.time(save(scRNA.norm.pca, file = "scRNA.norm.pca.Jack.Rdata"))


#load('./scRNA.norm.pca.Jack.Rdata')
## JackStrawPlot()函数：可视化比较每个主成分的 p 值分布和均匀分布（虚线）
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)

png("04.pca_cluster.png",w=10,h=7,units = 'in',res = 600,family='Times')
plot_pca
dev.off()
pdf("04.pca_cluster.pdf",w=10,h=7,family='Times')
plot_pca
dev.off()

plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)

pdf("05.pca_sd.pdf",w=6,h=5,family='Times')
plot_elbow
dev.off()

png("05.pca_sd.png",w=6,h=5,units = 'in',res = 600,family='Times')
plot_elbow
dev.off()

load('./scRNA.norm.pca.Jack.Rdata')
scRNA.norm.pca.clu <- FindNeighbors(scRNA.norm.pca, dims = 1:30)

## FindClusters()函数作用：以迭代方式将细胞分组在一起,resolution参数表示区分细胞群的分辨率，resolution越大，分的群也就越多。
scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu, resolution = 0.4) # 0.4


## UMAP聚类分析
UMAP <- RunUMAP(object = scRNA.norm.pca.clu, dims = 1:30)                      #UMAP聚类
system.time(save(UMAP, file = "UMAP_resolution_0.4.Rdata"))



load("UMAP_resolution_0.4.Rdata")
theme.set <-  theme(#axis.title.x=element_blank(),
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.position = 'right',legend.direction = 'vertical',
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))


pdf(file = paste0("06.UMAP1.pdf"),width = 8,height = 6,family="Times")
a <- dev.cur()
png(file = paste0("06.UMAP1.png"),width= 8, height= 6, units="in", res=600,family="Times")
dev.control("enable")
par(mar = c(2,2,2,2),cex=1.5,family="Times")
DimPlot(UMAP,reduction = 'umap',label = T)+
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +theme.set
dev.copy(which = a)
dev.off()
dev.off()

