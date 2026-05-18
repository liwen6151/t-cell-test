rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

options(stringsAsFactors = F) 
suppressPackageStartupMessages({
  library(rstudioapi)
  library(dplyr)
  library(GEOquery)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(data.table)
  library(stringr)
})
conflicted::conflict_prefer_matching("^filter$|^select$|^arrange$|^mutate$|", "dplyr",quiet = T)



# 04.cell annot  --------
print("# 04.cell annot-------")
dir <-file.path('04_cell_annot')
dir.create(dir,recursive = T,showWarnings = F)
setwd(dir)

load('../03_PCA_UMAP/UMAP_resolution_0.4.Rdata')

table(UMAP@meta.data$orig.ident)


scRNA <- UMAP



options(stringsAsFactors = F)
library(tidyverse)
library(patchwork)
theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))
if(T){
  ## 6.1 识别各celltype中的marker基因 ----
  #利用FindAllMarkers寻找每一种细胞类型的topMarker基因 
  # check the current active plan
  # library(future)
  # plan()
  # plan("multicore", workers = 128)
  # plan()
  print('start at:')
  print(Sys.time())
  DefaultAssay(scRNA) <- "RNA"
  scRNA <- JoinLayers(scRNA)
  all.markers <- FindAllMarkers(
    scRNA,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    test.use = 'wilcox',
    return.thresh = 0.01
  )
  print('FindAllMarkers:');print(Sys.time())
  write.csv(all.markers,file = '06.AllMarkers_cluster.csv',quote = F)
  saveRDS(all.markers,'all.markers_cluster.rds')
  #前台处理------
  all.markers <- readRDS('all.markers_cluster.rds')
  
  
  library(tidyverse)
  top <- all.markers %>% group_by(cluster) %>% top_n (n=3,wt=avg_log2FC) #提取每个cluster的top2基因
  write.csv(top,'06.topheat(cluster).csv')
  #load('../03_UMAP/scRNA_UMAP.Rdata')
  DefaultAssay(scRNA) <- "RNA"
  scRNA <- JoinLayers(scRNA)
  # 点阵图
  dotplot1 <- DotPlot(scRNA,features = unique(top$gene)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     color = "black",size=9))
  dotplot1
  ggsave(filename = '06.markers_dotplot(cluster).png',dotplot1,w=18,h=8,dpi = 600)
  ggsave(filename = '06.markers_dotplot(cluster).pdf',dotplot1,w=18,h=8)
  
  # 热图
  heatmap1 <- DoHeatmap(scRNA,features = top$gene,label = F)+
    labs(title="", x="Cells separated by clusters", y = "",size=40)
  
  ggsave(filename = '06.topheat(cluster).png',heatmap1,w=8,h=7)
  ggsave(filename = '06.topheat(cluster).pdf',heatmap1,w=8,h=7)
}



#文献注释
# 6.2 绘制文章注释细胞的Marker图 ----------------------------------------------------
## 01 查看marker gene 表达 ----
genes_to_check <- list(
  "Endothelial" = c("CLDN5", "SELE",'A2M') ,
  "Fibroblasts" = c("COL6A2", "DCN",'C11orf96'),
  "Epithelial" = c("WFDC2", "EPCAM",'CLDN3'),
  "Mast_Cell" = c("GATA2", "TPSAB1", "TPSB2"),
  "T_Cell" = c("CD3D","CD8A","CD3E", "GZMB","GNLY") ,
  "Macrophage" = c("HLA-DRB1", "FCER1G", "C15orf48")
)
unlist(genes_to_check)%>%.[duplicated(.)]
#0
genes<-unlist(genes_to_check)%>%unique()
genes[!genes%in%rownames(scRNA)] %>% print
#0
genes<-genes[genes%in%rownames(scRNA)]
p <- DotPlot(scRNA, features = genes_to_check,
             assay='RNA'  )   + theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                                                 color = "black",size=9))

p
ggsave('02.check_paper_markers.pdf',height = 9,width = 25)
ggsave('02.check_paper_markers.png',height = 9,width = 25)

genes<-genes[genes%in%rownames(scRNA)]
# X06 <- read.csv("06.topheat(cluster).csv")
# X06 %>% filter(cluster=='17') %>% pull(gene) %>% cat
# X06 %>% filter(cluster=='18') %>% pull(gene)%>% cat
# X06 %>% filter(cluster=='23') %>% pull(gene)%>% cat
# 创建注释列表
cluster_annotations <- list(
  "Endothelial" = c(0,17) %>% as.character(),
  "Fibroblasts" = c(1,2,3,7,8,11,14,20) %>% as.character(),
  "Epithelial" = c(5,10,12,13,18) %>% as.character(),
  "Mast_Cell" = c(16) %>% as.character(),
  "T_Cell" = c(4,6,15) %>% as.character(),
  "Macrophage" = c(9,19) %>% as.character()
) %>% unlist()

new.cluster.ids =gsub('\\d+$','',names(cluster_annotations))
names(new.cluster.ids)=cluster_annotations
UMAP<- RenameIdents(scRNA, new.cluster.ids)
UMAP$celltype = Idents(UMAP)

saveRDS(UMAP,'UMAP.celltype.rds')

table(UMAP$celltype) %>% sort() %>% print


# 定义完整的颜色 palette（高辨识度版，覆盖所有细胞类型）
global_palette <- c(
  "Endothelial" = "#377EB8",         # 蓝色
  "Fibroblasts" = "#4DAF4A",         # 绿色
  "Epithelial" = "#984EA3",    # 紫色
  "Mast_Cell" = "#E41A1C",        # 红色
  "T_Cell" = "#FF7F00",             # 橙色
  "Macrophage" = "#FFEAA7"       # 灰色
)



library(RColorBrewer)
library(ggsci)
library(grDevices)
library(IOBR)
p1 <- DotPlot(UMAP, features = genes,group.by = 'celltype',
              assay='RNA'  )   + theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                                                  color = "black",size=9))

ggsave('02.check_paper_markers_dot.pdf',p1,height = 5,width = 12)
ggsave('02.check_paper_markers_dot.png',p1,height = 5,width = 12)
anno_plot2 <-DimPlot(UMAP, reduction = "umap", group.by = "celltype", 
                     label = F,label.size = 6,
                     cols = global_palette)+theme.set

anno_plot2
ggsave(filename = '02.check_paper_markers_celltype.pdf',anno_plot2,w=10,h=6)
ggsave(filename = '02.check_paper_markers_celltype.png',anno_plot2,w=10,h=6,dpi = 600)

##4.2.可视化细胞的上述比例情况（根据细胞类型及进行可视化） -----------------------------
#去掉小提琴图上不合适的基因
#color2=palettes(category = "random",length(genes),show_col=F)
pv<-VlnPlot(UMAP,
            features = genes,
            split.by = 'celltype',stack = TRUE,
            flip=F,
            cols=global_palette)+
  theme.set+theme(axis.ticks.x = element_blank(), 
                  axis.text.x = element_blank(),
                  plot.margin = ggplot2::margin(t=0.1,r=0.1,b=0.1,l=0.5,unit = 'in')
  )
pv
ggsave(filename = '06.2.check_paper_markers_exp.pdf',pv,w=12,h=6)
ggsave(filename = '06.2.check_paper_markers_exp.png',pv,w=12,h=6,dpi = 300)
#celltype umap-------
library(IOBR)
anno_plot3 <-DimPlot(UMAP, reduction = "umap", group.by = "celltype", #split.by = 'group',
                     label = T,label.size = 6,repel = T,
                     cols = global_palette)+theme.set+
  theme(text = element_text(family = "Times",size = 10))

anno_plot3
ggsave("05.3.umap.celltype.pdf", plot = anno_plot3, width = 10, height = 8)
ggsave("05.3.umap.celltype.png", plot = anno_plot3, width = 12, height = 8,dpi = 300)


