rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("02_VariableFeatures")) {dir.create("02_VariableFeatures")}
setwd("02_VariableFeatures/")


# Seurat v5的标准整合流程-----------------------------
load('../01_QC/scRNA_qc.Rdata')

library(Seurat)
library(tidyverse)

table(scRNA$orig.ident)

scRNA <- JoinLayers(scRNA)
###########表达量标准化
scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##########寻找高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(scRNA), 10)
# plot variable features with and without labels
p <- VariableFeaturePlot(scRNA)
p <- LabelPoints(plot = p, 
                 points = top10, 
                 repel = T)
png("01.feature_selection.png",w=9,h=6,units = 'in',res = 600,family='Times')
p
dev.off()
pdf("01.feature_selection.pdf",w=9,h=6,family='Times')
p
dev.off()


#####表达量scale处理
scRNA <- ScaleData(scRNA)


system.time(save(scRNA, file = "scRNA_VariableFeatures.Rdata"))

