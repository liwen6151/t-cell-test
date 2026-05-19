rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (!dir.exists("08_GSEA")){dir.create("08_GSEA")}
setwd("08_GSEA")

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(DESeq2)
library(tidyverse)
library(lance)
## 
dat <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', check.names = F, row.names = 1)%>%lc.tableToNum()
colnames(dat) <- gsub('.', '-', colnames(dat), fixed = T)
dat <- round(dat, digits = 0) # round是舍入函数 


# 读取gmt文件中的基因列表
gmt <- read.gmt("/data/nas1/zhangzhaolei/pipline/GSEA/c2.cp.kegg.v7.4.symbols.gmt") #读gmt文件
length(unique(gmt$term))

# 用substr函数在TCGA数据中提取样本信息
group <- read.csv("../05_Prognostic_model/train/risk.csv")

names(group)[names(group) == "X"] <- "Sample" # 重命名列Sample为Gene_Symbol
group$risk <- factor(ifelse(group$risk == "high risk", "High", "Low"))
levels(group) = c('High','Low')
table(group$risk)

names(group)[names(group) == "risk"] <- "Group" # 重命名列risk为Group
table(group)
data1 <- dat[, group$Sample]

identical(colnames(data1),group$Sample)


dds <- DESeqDataSetFromMatrix(countData = data1, colData = group, design = ~Group)
dds <- dds[rownames(counts(dds)) > 1, ] 
dds <- estimateSizeFactors(dds) 
##提取标准化后的数据 
normalized_counts <- counts(dds, normalized = T)
dds <- DESeq(dds)
## 提取差异结果
res <- results(dds, contrast = c("Group","High", "Low"))
res <- res[order(res$padj), ]
head(res)
summary(res)

DEG <- as.data.frame(res)
DEG <- na.omit(DEG)

DEG_write <- cbind(GeneSymbol = rownames(DEG), DEG)
write.csv(DEG_write, file = "TCGA_DEGs_risk_all_DEGs.csv",quote = F,row.names = F)


# GSEA分析----------------------------------------------
tempOutput <- read.csv('TCGA_DEGs_risk_all_DEGs.csv', row.names=1)
head(tempOutput, n=3)
genelist <- data.frame(Gene_Symbol = rownames(tempOutput), logFC = tempOutput$log2FoldChange)
#开始ID转换
gene <- bitr(genelist$Gene_Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(ENTREZID=gene$ENTREZID ,#可以是foldchange
                      Gene_Symbol = gene$SYMBOL) #记住你的基因表头名字
gene_df <- merge(gene_df,genelist,by="Gene_Symbol")
geneList <- gene_df$logFC #第二列可以是folodchange，也可以是logFC
names(geneList) <- gene_df$Gene_Symbol #使用转换好的ID
geneList <- sort(geneList, decreasing = T) #从高到低排序
set.seed(1)

## GSEA分析
res_GSEA <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 0.05, eps = 0)
sortGESA <- data.frame(res_GSEA)
sortGESA <- sortGESA[order(sortGESA$p.adjust, decreasing = F),]#按照enrichment score从高到低排序
if(length(sortGESA$ID) > 5) {
  paths <- rownames(sortGESA[c(1: 5), ]) # 选取需要展示的通路ID（比如前5条）
}else {
  paths <- sortGESA$ID
}
## 保存GSEA结果
write.csv(sortGESA, file = "TCGA_gsea_res.csv")

## 绘图
library(GseaVis)
library(psych)
library(reshape2)
p <- gseaNb(object = res_GSEA,
            geneSetID = paths,
            subPlot = 2,
            termWidth = 45,
            # legend.position = c(0.72,0.8),
            addPval = F,
            rmHt = F,
            # pvalX = 0.99,
            # pvalY = 0.99, 
            newGsea = F,
            curveCol = c('#FFDAB9', '#00FFFF', '#F8766D', '#959897', '#90EE90','#00A9FF','darkgoldenrod1','#E01516','#35A132','black'))
p
pdf(file = '01.Multi_paths_gsea.pdf', family = "Times", w=10, h=7)
print(p)
dev.off()
png(file = '01.Multi_paths_gsea.png', family = "Times", w=10, h=7, units = 'in', res = 600)
print(p)
dev.off()

