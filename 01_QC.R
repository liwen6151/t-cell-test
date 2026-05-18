rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")


if (!dir.exists("01_QC")) {dir.create("01_QC")}
setwd("01_QC/")


load('../00_scRNA_rawdata/scRNA_orig.Rdata')



## 01 QC 1 可视化质控前的质控指标情况-----
feats <- c("nFeature_RNA", "nCount_RNA","percent.mt") #
p1=VlnPlot(scRNA, group.by = "orig.ident",add.noise = F, features = feats, pt.size = 0, ncol = 3) + 
  NoLegend() #绘制小提琴图，指定feats； pt.size代表点的大小；ncol代表2列
p1
library(ggplot2)
ggsave("01.vlnplot_before_qc.pdf", plot = p1, width = 12, height = 5)
ggsave("01.vlnplot_before_qc.png", plot = p1, width = 12, height = 5,dpi = 600)
write.csv(list(cell=length(colnames(scRNA)),gene=length(rownames(scRNA))),'01.before_qc_data.csv')

## 绘制nFeatrey与nCount两个指标之间的相关性图；越高代表测序越好
p3=FeatureScatter(scRNA, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3,w=8,h=5)


## 01 QC 2 根据上述5个质控指标，过滤低质量细胞/基因-------
### 质控后 ###
### 设置质控标准
minGene=quantile(scRNA$nFeature_RNA,.02)
maxGene=quantile(scRNA$nFeature_RNA,0.95)
maxUMI=quantile(scRNA$nCount_RNA,.95)
minUMI=quantile(scRNA$nCount_RNA,.02)
minGene%>% print;maxGene%>% print;maxUMI%>% print;minUMI%>% print

scRNA.org<-scRNA
minGene=200
maxGene=5000 #基因数
maxUMI=25000 #表达量之和分子数
minUMI=200

#length(colnames(scRNA))   #cell
#length(rownames(scRNA)) #gene
dim(scRNA.org)%>% print
#[1] 26145 35868
pctMT=20
scRNA <- subset(scRNA.org, subset =
                  nCount_RNA < maxUMI &
                  nCount_RNA > minUMI & 
                  nFeature_RNA > minGene &
                  nFeature_RNA < maxGene &
                  percent.mt < pctMT )
dim(scRNA)%>% print
#[1] 26145 31303
# 质控后小提琴图
#feats <- c("nFeature_RNA", "nCount_RNA") #, "percent.mt"没有线粒体基因
p2=VlnPlot(scRNA, group.by = "orig.ident",add.noise = F, features = feats, pt.size = 0, ncol = 3) + 
  NoLegend() #绘制小提琴图，指定feats； pt.size代表点的大小；ncol代表2列
p2
library(ggplot2)
ggsave("02.vlnplot_after_qc.pdf", plot = p2, width = 12, height = 5)
ggsave("02.vlnplot_after_qc.png", plot = p2, width = 12, height = 5,dpi = 300)
write.csv(list(cell=length(colnames(scRNA)),gene=length(rownames(scRNA))),'02.after_qc_data.csv')

system.time(save(scRNA, file = "scRNA_qc.Rdata"))
