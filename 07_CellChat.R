
###细胞通讯
# # 细胞通讯1:--------
rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")


if (! dir.exists("./07_CellChat")){
  dir.create("./07_CellChat")
}
setwd("./07_CellChat")


library(CellChat)
library(Seurat)
library(ggplot2)
# #读入seurat处理后的rds文件
scRNA <-readRDS('../04_cell_annot/UMAP.celltype.rds')


data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA$data # Controlized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data

meta$group<-"UCEC"
table(meta$group)

# ##疾病和对照分开画（疾病组UCEC）
cell.use = rownames(meta)[meta$group == "UCEC"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)
table(meta$orig.ident)

# ##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

# ##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <-CellChatDB
# ##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# # project gene expression data onto PPI network (optional)
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <-projectData(cellchat, PPI.human)

# ##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
# ##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)
df.net<-na.omit(df.net)
write.csv(df.net,'01.df.net.UCEC.csv',row.names = F,quote = F)
# ##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# ##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

# #### circle
pdf('01.net_number.UCEC.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions-UCEC")
dev.off()

png('01.net_number.UCEC.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions-UCEC")
dev.off()

pdf('02.net_weight.UCEC.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weight/strength-UCEC")
dev.off()

png('02.net_weight.UCEC.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weight/strength-UCEC")
dev.off()
levels(cellchat@idents)


#######单个细胞与其它细胞的相互作用强度
mat <- cellchat@net$weight
pdf('3-1.single_circle_weight.UCEC.pdf', w=15, h=10)
par(mfrow=c(1, 1), xpd=T)  # 每行放4个图
par(mar=c(3, 3, 3, 3))      # 设置边距
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

png('3-1.single_circle_weight.UCEC.png',w=15,h=10,units='in',res=600)
par(mfrow=c(3, 3), xpd=T)  # 每行放4个图
par(mar=c(3, 3, 3, 3))      # 设置边距
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#######单个细胞与其它细胞的相互作用数量
mat <- cellchat@net$count
pdf('3-2.single_circle_number.UCEC.pdf', w=15, h=10)
par(mfrow=c(1, 1), xpd=T)  # 每行放4个图
par(mar=c(3, 3, 3, 3))      # 设置边距
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

png('3-2.single_circle_number.UCEC.png',w=15,h=10,units='in',res=600)
par(mfrow=c(3, 3), xpd=T)  # 每行放4个图
par(mar=c(3, 3, 3, 3))      # 设置边距
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat<- netAnalysis_computeCentrality(cellchat)

##关键细胞受体配体信号气泡图--Macrophage
table(cellchat@meta$celltype)
a<- as.data.frame(cellchat@var.features[["features.info"]])
range(a$pvalues)
range(a$logFC)
library(openxlsx)
write.xlsx(a, file = "04.celltype_features(p0.05).xlsx")

p1 <- netVisual_bubble(cellchat,
                       sources.use = c("Epithelial"),
                       targets.use = c(  "Endothelial",
                                         "Fibroblasts",     
                                         "Epithelial",  
                                         "Mast_Cell",      
                                         "T_Cell",          
                                         "Macrophage"),
                       angle.x = 45,
                       remove.isolate=FALSE)
##调整 Y 轴字体大小
p1 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

##调整 Y 轴字体大小
p1 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

ggsave(filename="04.Epithelial_UCEC.pdf",p1,height=10,width=8,dpi=600)
ggsave(filename="04.Epithelial_UCEC.png",p1,height=10,width=8,dpi=600)




p2 <- netVisual_bubble(cellchat,
                       targets.use = c("Epithelial"),
                       sources.use = c(  "Endothelial",
                                         "Fibroblasts",     
                                         "Epithelial",  
                                         "Mast_Cell",      
                                         "T_Cell",          
                                         "Macrophage"),
                       angle.x = 45,
                       remove.isolate=FALSE)
##调整 Y 轴字体大小
p2 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

##调整 Y 轴字体大小
p2 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

ggsave(filename="05.Epithelial_UCEC.pdf",p2,height=10,width=8,dpi=600)
ggsave(filename="05.Epithelial_UCEC.png",p2,height=10,width=8,dpi=600)
save(cellchat,file = 'UCEC_cellchat.RData')



p3 <- netVisual_bubble(cellchat,
                       sources.use = c("Epithelial"),
                       targets.use = c( "T_Cell"),
                       angle.x = 45,
                       remove.isolate=FALSE)
##调整 Y 轴字体大小
p3 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

##调整 Y 轴字体大小
p3 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))



p4 <- netVisual_bubble(cellchat,
                       targets.use = c("Epithelial"),
                       sources.use = c( "T_Cell"),
                       angle.x = 45,
                       remove.isolate=FALSE)
##调整 Y 轴字体大小
p4 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

##调整 Y 轴字体大小
p4 + theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))

plotc <- p3|p4
plotc

ggsave(filename="06.Epithelial_UCEC.pdf",plotc,height=8,width=8,dpi=600)
ggsave(filename="06.Epithelial_UCEC.png",plotc,height=8,width=8,dpi=600)
