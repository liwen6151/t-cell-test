# 基因突变景观瀑布图---------------------------
rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if(!dir.exists('./06_mutate')){
  dir.create('./06_mutate')
}
setwd('./06_mutate/')

library(TCGAmutations)
library(tidyverse)

maf1 <- TCGAmutations::tcga_load(study = "UCEC")

maf <- maf1
all.maf <- maf1

## 计算tmb值
x = tmb(maf = maf1)

library(maftools)
riskscore <- read.csv("../05_Prognostic_model/train/risk.csv", header = T, row.names = 1)
riskscore$sample <- rownames(riskscore)
riskscore$sample <- gsub('.', '-', riskscore$sample, fixed = TRUE)
High.sample <- riskscore$sample[which(riskscore$risk == "low risk")]
Low.sample <- riskscore$sample[which(riskscore$risk == "high risk")]


maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)


maf_sample  <- maf_sample[maf_sample$sample %in% riskscore$sample,]
maf_sample <- merge(maf_sample, riskscore, by = 'sample')

sample <- subset(maf_sample)$barcode
maf <- subsetMaf(maf, tsb = sample)


##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample, ]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "01.oncoplot.high.pdf", height = 8, width = 10)
oncoplot(maf = maf.high, top = 20)
dev.off()

png(file = "01.oncoplot.high.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.high, top = 20)
dev.off()


low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "02.oncoplot.low.pdf", height = 8, width = 10)
oncoplot(maf = maf.low, top = 20)
dev.off()

png(file = "02.oncoplot.low.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.low, top = 20)
dev.off()

