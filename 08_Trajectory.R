
###########提取关键细胞重聚类达
rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("08_Trajectory")) {dir.create("08_Trajectory")}
setwd("08_Trajectory/")


library(GEOquery)
library(tidyverse)
library(CellChat)
library(dplyr)
library(Seurat)
library(hdf5r)
library(ggpubr)
library(data.table)
library(CellChat)
library(Seurat)
library(ggplot2)

UMAP <- readRDS('../04_cell_annot/UMAP.celltype.rds')


table(UMAP@meta.data$celltype)


##########提取细胞进行重聚类

scRNA <- UMAP[, UMAP@meta.data$celltype %in% c('Epithelial')]
scRNA <- JoinLayers(scRNA)   #需要先将数据合并
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)

FindVariableFeatures(scRNA,
                     selection.method = 'vst',
                     nfeatures = 2000)
scRNA <- ScaleData(scRNA)
scRNA <-
  RunPCA(scRNA,
         features = VariableFeatures(object = scRNA))

scRNA <-
  JackStraw(scRNA,
            num.replicate = 100,
            dims = 50)
scRNA <-
  ScoreJackStraw(scRNA, dims = 1:50)

p1 <- JackStrawPlot(scRNA, dims = 1:50)
ggsave('01.JackStrawPlot.png',
       p1,
       w = 10,
       h = 5)
ggsave('01.JackStrawPlot.pdf',
       p1,
       w = 10,
       h = 5)

p2 <- ElbowPlot(scRNA, ndims = 50)
ggsave('02.subtype_ElbowPlot.png', p2, w = 6, h = 5)
ggsave('02.subtype_ElbowPlot.pdf', p2, w = 6, h = 5)


scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.1)

######按照样本分开展示umap图
colors <-c("#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9",
           "#ebce8e","#f9766e","#4e79a6","#7587b1","#c7deef",
           "#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
           "#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
           "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
           "#279772","#add387","#d9b71a","#fbbab6","#e97371",
           "#e1c548","#f0e2a3","#aedd2f","#d7ee96","#a199be",
           "#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#e5c06e",
           "#bac4d0","#45337f")

p3 <-
  DimPlot(
    scRNA,
    reduction = "umap",
    cols = colors,
    group.by = "seurat_clusters",
    label = T
  ) +
  ggtitle('Epithelial')
p3
ggsave(filename = '03.subtype_umap.png', p3, w = 6, h = 5)
ggsave(filename = '03.subtype_umap.pdf', p3, w = 6, h = 5)

system.time(save(scRNA, file = "Epithelial.Rdata"))  


#####关键细胞拟时序分析-- Memory CD8+ T cells
rm(list = ls())

library(monocle)
library(Seurat)
library(tidyverse)

#trace("project2MST", where = asNamespace("monocle"), edit = TRUE)  #运行一次即可
#将if(class(projection) != 'matrix')projection <- as.matrix(projection)修改为projection <- as.matrix(projection)

load('Epithelial.Rdata') #挑最显著细胞簇!!!
table(scRNA$seurat_clusters)
Mono_tj = scRNA
table(scRNA$seurat_clusters)

Mono_matrix = GetAssayData(Mono_tj, layer = "count", assay = "RNA")
feature_ann = data.frame(gene_id = rownames(Mono_matrix),
                         gene_short_name = rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Cluster"
Mono_pd = new("AnnotatedDataFrame", data = sample_ann)
Mono.cds = newCellDataSet(
  Mono_matrix,
  phenoData = Mono_pd,
  featureData = Mono_fd,
  expressionFamily = negbinomial.size()
)
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()

#plot_ordering_genes(Mono.cds)
Mono.cds = reduceDimension(
  Mono.cds,
  max_components = 2,
  verbose = T,
  norm_method = "log"#, residualModelFormulaStr = "~nFeature_RNA+group"
)
save(Mono.cds, file = "cds.RData")

####orderCells报错if (class(projection) != "matrix") projection <- as.matrix(projection) : the condition has length > 1
##可能是monocle包不行，我这里选择本地安装2.24版本的monocle包后成功运行了
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "orderCells.cds.RData")
Mono.cds$Subtype<- Mono.cds$seurat_clusters
p = plot_cell_trajectory(
  Mono.cds,
  color_by = "Pseudotime",
  cell_size = 1,
  theta = 180,
  size = 1,
  show_backbone = TRUE,
  show_branch_points = F
) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
ggsave(
  "04.Pseudotime_cell.png",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)
ggsave(
  "04.Pseudotime_cell.pdf",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

p3 = plot_cell_trajectory(
  Mono.cds,
  color_by = "State",
  cell_size = 1,
  theta = 180,
  size = 1,
  show_backbone = TRUE,
  show_branch_points = F
) +
  scale_color_manual(
    values = c(
      "#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9",
      "#ebce8e","#f9766e","#4e79a6","#7587b1","#c7deef",
      "#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
      "#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
      "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
      "#279772","#add387","#d9b71a","#fbbab6","#e97371",
      "#e1c548","#f0e2a3","#aedd2f","#d7ee96","#a199be",
      "#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#e5c06e",
      "#bac4d0","#45337f"
    )
  ) +
  theme(text = element_text(size = 14))
p3
ggsave(
  "05_State_Trajectory.png",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)
ggsave(
  "05_State_Trajectory.pdf",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

p5 = plot_cell_trajectory(
  Mono.cds,
  color_by = "Subtype",
  cell_size = 1,
  theta = 180,
  size = 1,
  show_backbone = TRUE,
  show_branch_points = F
) +
  scale_color_manual(
    values = c(
      "#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9",
      "#ebce8e","#f9766e","#4e79a6","#7587b1","#c7deef",
      "#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
      "#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
      "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
      "#279772","#add387","#d9b71a","#fbbab6","#e97371",
      "#e1c548","#f0e2a3","#aedd2f","#d7ee96","#a199be",
      "#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#e5c06e",
      "#bac4d0","#45337f"
    )
  ) +
  theme(text = element_text(size = 14))
p5
ggsave(
  filename = '06_Subtype_Trajectory.pdf',
  p5,
  w = 6,
  h = 6,
  units = "in",
  dpi = 300
)
ggsave(
  filename = '06_Subtype_Trajectory.png',
  p5,
  w = 6,
  h = 6,
  units = "in",
  dpi = 300
)

####关键基因表达

##hub基因表达-----



#########组合图
colors <-c("#1F77B4","#FF7F0E","#2CA02C","#3C5488FF","#925E9FFF","#FED439FF","#8491B4FF","#7E6148FF","#FDAF91FF","#B24745FF","#6699FFFF","#99991EFF","#FFCCCCFF",
           "#99CCFFFF","#FFCC00FF","#8C564BFF","#BCBD22FF","#996600FF","#5CB85CFF","#F39B7FFF","#CE3D32FF","#749B58FF",
           "#466983FF","#F0E685FF","#D595A7FF","#924822FF","#7A65A5FF","#C75127FF","#FFA319FF","#8A9045FF","#8F3931FF","#00AF66FF",
           "#748AA6FF","#D0DFE6FF","#C71000FF","#008EA0FF","#8A4198FF","#D5E4A2FF","#5A9599FF","#FF6348FF","#B7E4F9FF","#FF95A8FF",
           "#526E2DFF","#FB6467FF","#E89242FF","#69C8ECFF","#917C5DFF","#709AE1FF","#D2AF81FF","#FD7446FF","#98df8aFF")


s.genes <-  c("SPDEF","CENPF","CDKN2A","E2F1","KLF2","TUBB2B","CKMT1B")#"OTC"
s.genes %in% rownames(Mono.cds)
#cds_subset <- Mono.cds[s.genes,]


p <-plot_pseudotime_heatmap(Mono.cds[s.genes,],
                            num_clusters = 4, 
                            show_rownames = TRUE,
                            cores = 4,return_heatmap = TRUE,
                            hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
p

ggsave('07.热图_gene.png',p,w=4,h=3, units = "in", dpi = 300)
ggsave('07.热图_gene.pdf',p,w=4,h=3, units = "in", dpi = 300)


save.image("hubgene_trajectory.RData")

