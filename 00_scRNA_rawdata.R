rm(list = ls())
setwd("/data/nas1/lijia/124_KYGW-61101-6-NKT-KM114/")

if (!dir.exists("00_scRNA_rawdata")) {dir.create("00_scRNA_rawdata")}
setwd("00_scRNA_rawdata/")


library(Seurat)
library(dplyr)
library(Matrix)
library(scales)
library(harmony)
library(ggplot2)
#library(GSEABase)
library(tidyr)
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(data.table)


#scRNA数据预处理---------------

fs=list.files('../00_scRNA_rawdata/GSE173682_RAW/')
dir <-paste0("GSE173682_RAW/",fs)
dir
samples_name <- c("GSM5276933","GSM5276934","GSM5276935","GSM5276936","GSM5276937")

scRNAlist <- list() 
for(i in 1:length(dir)) {
  counts <- Read10X(dir[i]) 
  #不设置min.cells过滤基因会导致CellCycleScoring报错：
  scRNAlist[[i]] <-CreateSeuratObject(
    counts,
    project = samples_name[i],
    min.cells = 3,
    min.features = 200
  )
  scRNAlist[[i]]@meta.data$orig.ident <- samples_name[i]
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <-RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
  #计算线粒体基因比例
  if (T) {
    scRNAlist[[i]][["percent.mt"]] <-PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  }
  #计算核糖体基因比例
  if (T) {
    scRNAlist[[i]][["percent.rb"]] <-PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #计算红细胞基因比例
  # if(T){ HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  # HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
  # scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) }
}

names(scRNAlist) <- samples_name 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)]) 
scRNA


table(scRNA@meta.data$orig.ident)


save(scRNA,file = 'scRNA_orig.Rdata')


