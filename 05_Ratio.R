rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("05_Ratio")) {dir.create("05_Ratio")}
setwd("05_Ratio/")


####标记基因在差异细胞中的表达------------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(ggsci)
library(paletteer)
library(patchwork)
library(ComplexHeatmap)
library(circlize)



UMAP<-readRDS('../04_cell_annot/UMAP.celltype.rds')

table(UMAP$celltype)

gene <-  c("SPDEF","CENPF","CDKN2A","E2F1","KLF2","TUBB2B","CKMT1B","OTC")


library(RColorBrewer)

pal <- rev(brewer.pal(n = 10, name = "RdBu"))
gene_colors <- pal[1:4]


PMgenes <- gene

pdf(file = paste0( "03.hub_gene_point_map.pdf"), width = 7, height = 5)
print(DotPlot(UMAP, features = PMgenes, group.by = "celltype",assay='RNA')+
        theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle = 45))+
        labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) + scale_size(range =c (3, 7)) +
        scale_color_gradientn(values = seq(0,1,0.1),colours = c('#330066','#336699','#66CC66','#FFCC33')))
dev.off()
# 
# 
png(file = paste0( "03.hub_gene_point_map.png"), width = 7, height = 5,res = 600,units = 'in')
print(DotPlot(UMAP, features = PMgenes, group.by = "celltype",assay='RNA')+
        theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle = 45))+
        labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+ scale_size(range =c (3, 7)) +
        scale_color_gradientn(values = seq(0,1,0.1),colours = c('#330066','#336699','#66CC66','#FFCC33')))
dev.off()


