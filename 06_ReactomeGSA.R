rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("06_ReactomeGSA")) {dir.create("06_ReactomeGSA")}
setwd("06_ReactomeGSA/")


options(stringsAsFactors = F)

library(pheatmap)
library(CellChat) 
library(Seurat)
library(ggplot2)
library(ReactomeGSA)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)

UMAP<-readRDS('../04_cell_annot/UMAP.celltype.rds')

table(UMAP$celltype)



##添加细胞类型信息
pbmc <- UMAP


##查看meta信息
head(pbmc@meta.data)



##将细胞类型信息赋值给我们的Seurat对象pbmc
Idents(pbmc) <- "celltype"
##调用analyse_sc_clusters函数，该函数对scRNA对象进行GSVA分析，并返回一个gsva_result对象，该对象是一个GSVAList对象，包含了分析的结果。同时打印了一些信息，显示分析的进度。
options(reactome_gsa.url = "http://gsa.reactome.org")
gsva_result <- analyse_sc_clusters(pbmc, verbose = TRUE)
save(gsva_result,file='gsva_result.RData')


##从gsva_result对象中提取了基因集在不同簇之间的倍数变化，并显示了数据框的前六行。
head(gsva_result@results$Seurat$fold_changes)
##从gsva_result对象中提取了基因集表达矩阵，并去掉了列名中的".Seurat"后缀。只显示矩阵的前三行。
pathway_expression <- pathways(gsva_result)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
pathway_expression[1:3,]
##找到最大差异表达的基因集，通过计算每个基因集的最大和最小表达值之间的差值。
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
max_difference$diff <- max_difference$max - max_difference$min
##我们按照差异的降序对基因集进行排列，并显示数据框的前15行。
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
head(max_difference,15)


plot_num = 20
plot_gsva <- pathway_expression[rownames(max_difference[1:plot_num,]),]
write.csv(plot_gsva,file = "gsva.csv")
pheatmap(
  filename = "INS.max_difference_top20.png",
  #height = 9,width = 12,
  plot_gsva[,-1],   # 使用数据矩阵，不转置
  scale = "column",       # 对列进行标准化，因为现在列是原始的行
  angle_col = 45,      # 设置列名角度为45度
  cellwidth = 60,      # 设置每个单元格的宽度
  cellheight = 30,     # 设置每个单元格的高度
  labels_row = plot_gsva[,1],  # 使用数据的第一列作为行标签
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)  # 定义颜色渐变
)

pheatmap(
  filename = "INS.max_difference_top20.pdf",
  #height = 9,width = 12,
  plot_gsva[,-1],   # 使用数据矩阵，不转置
  scale = "column",       # 对列进行标准化，因为现在列是原始的行
  angle_col = 45,      # 设置列名角度为45度
  cellwidth = 60,      # 设置每个单元格的宽度
  cellheight = 30,     # 设置每个单元格的高度
  labels_row = plot_gsva[,1],  # 使用数据的第一列作为行标签
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)  # 定义颜色渐变
)


