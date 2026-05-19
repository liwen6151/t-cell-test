##tcga-ucec-----------------
rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if(!dir.exists('./00_rawdata')){
  dir.create('./00_rawdata')
}
setwd('./00_rawdata/')

if (!dir.exists("TCGA")) {dir.create("TCGA")}
setwd("TCGA")

# ###TCGA--------
library(readr)
library(readxl)
library(tidyverse)
library(tibble)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)

## 下载日期：2025.7.16
# query <- GDCquery(
#    project = "TCGA-UCEC",
#    data.category = "Transcriptome Profiling",
#    data.type = "Gene Expression Quantification",
#    workflow.type = "STAR - Counts"
# )
# length(query$results[[1]]$file_id)
# GDCdownload(query, directory = "GDCdata", method = "api", files.per.chunk = 20)
# expquery <- GDCprepare(query,directory = "GDCdata",summarizedExperiment = T)
# save(expquery,file = "UCEC.gdc_2025.7.rda")
# 
# 
# query <- GDCquery(
#    project = "TCGA-UCEC",
#    data.category = "Clinical",
#    data.type = "Clinical Supplement",
#    data.format = "BCR XML"
# 
# )
# GDCdownload(query, directory = "GDCdata", method = "api", files.per.chunk = 20)
# clinical <- GDCprepare_clinic(query, clinical.info = "patient")
# save(clinical,file = "UCEC.gdc.clinical_2025.7.rda")
# print('finished')

load('./UCEC.gdc_2025.7.rda')
rowdata <- rowData(expquery) %>% as.data.frame()
table(rowdata$gene_type)

probe2symbol <- rowdata %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  dplyr::select(ID = gene_id, symbol = gene_name)
sample <- data.frame(sample = expquery$sample,
                     barcode = expquery$barcode)
duplicated_samples1 <- sample %>%
  group_by(sample) %>%
  filter(n() == 1)
duplicated_samples2 <- sample %>%
  group_by(sample) %>%
  filter(n() > 1)
if(nrow(duplicated_samples2)>1){
  duplicated_samples2 <- duplicated_samples2%>%
    arrange(sample,
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,5]),
            desc(str_split(duplicated_samples2$barcode, "-", simplify = TRUE)[,6]))
  duplicated_samples2 <- duplicated_samples2[!duplicated(duplicated_samples2$sample),]
}
barcode <- c(duplicated_samples1$barcode,duplicated_samples2$barcode)


## 删除重复样本
count <- setNames(as.data.frame(assay(expquery, "unstranded")), expquery$barcode)
count <- count[,barcode]
colnames(count) <- substr(colnames(count), 1, 16)

#count <- count[,!duplicated(colnames(count))]
range(count)
dat_count <- count
dat_count$ID <- rownames(dat_count) %>% as.character()

##id转换
dat_count <- dat_count %>%
  inner_join(probe2symbol, by = 'ID') %>%
  dplyr::select(ID, symbol, everything()) %>% 
  mutate(rowMean = rowMeans(.[grep('TCGA', names(.))])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol, .keep_all = T) %>% 
  dplyr::select(-rowMean)  

probeid <- dat_count$ID
rownames(dat_count) <- dat_count$symbol
dat_count[,c('ID', 'symbol')] <- NULL
dim(dat_count) 
#[1] 19938   585


## 样本分组
phenotype <- colData(expquery) %>% as.data.frame()
group <- data.frame(sample = phenotype$sample,
                    sample_type = phenotype$sample_type)
group <- group[!duplicated(group),]
group <- group[substr(group$sample, 16, 16) == "A", ]
table(group$sample_type)
group <- group %>%
  subset(sample_type %in% c('Primary Tumor', 'Solid Tissue Normal')) %>% 
  transform(Type = ifelse(sample_type == 'Solid Tissue Normal', 'Normal', 'Tumor')) %>%
  with(data.frame(row.names = sample, Type))
table(group$Type)
#Normal  Tumor 
# 35    539 
sample <- Reduce(intersect, list(rownames(group), colnames(dat_count)))
group.write0 <- group[sample,,drop=F] %>% dplyr::arrange(Type)
table(group.write0$Type)
# Normal  Tumor 
# 35    539 


disease <- 'UCEC'

write.csv(group.write0, file = paste0('./01.group_TCGA.',disease,'.csv'), quote=F) 

count.write <- dat_count[rownames(group.write0)]
identical(colnames(count.write), rownames(group.write0))
dim(count.write)
#19938   574
write.csv(count.write,file = paste0('./01.count_TCGA.',disease,'_mRNA.csv'), quote = F, row.names = T)


# ## 保留有生存数据的样本
survival <- subset(phenotype,select =c(sample,vital_status,days_to_death,days_to_last_follow_up))
survival <- survival[!duplicated(survival),]
table(survival$vital_status)
#survival <- survival[survival$vital_status != 'Not Reported',]
survival <- survival %>%
  transmute(sample,
            OS = as.integer(vital_status == "Dead"),
            OS.time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)) %>%
  filter(!is.na(OS.time) & OS.time > 0)
survival <- data.frame(row.names = survival$sample,
                       OS = survival$OS,
                       OS.time = survival$OS.time)

tumor <- rownames(group.write0[group.write0$Type == 'Tumor', , drop = FALSE])
survival.sample <- Reduce(intersect, list(tumor, rownames(survival)))
survival.write <- survival[survival.sample,]
dim(survival.write)
#491      2
write.csv(survival.write, file = paste0('./01.survival_TCGA.',disease,'.csv'), quote=F) 



##fpkm
fpkm <- setNames(as.data.frame(assay(expquery, "fpkm_unstrand")), expquery$barcode)
fpkm <- fpkm[probeid,barcode]
colnames(fpkm) <- substr(colnames(fpkm), 1, 16)
#fpkm <- fpkm[,!duplicated(colnames(fpkm))]
range(fpkm)
dat_fpkm <- log2(fpkm + 1)
range(dat_fpkm)
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm <- dat_fpkm %>%
  inner_join(probe2symbol, by = 'ID')
rownames(dat_fpkm) <- dat_fpkm$symbol
dat_fpkm[,c('ID', 'symbol')] <- NULL
dim(dat_fpkm) 
## 19938   585

fpkm.write <- dat_fpkm[rownames(count.write), colnames(count.write)]
identical(colnames(fpkm.write), rownames(group.write0))
identical(colnames(fpkm.write), colnames(count.write))
identical(rownames(fpkm.write), rownames(count.write))
dim(fpkm.write)
#  19938   574
write.csv(fpkm.write, file = paste0('./01.fpkmlog2_TCGA.',disease,'_mRNA.csv'), row.names = T, quote = F)

load('./UCEC.gdc.clinical_2025.7.rda')
write.csv(clinical, './data_clinical.csv') 



# 2. GSE119041(芯片数据)(验证集)---------------------------------------
rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if(!dir.exists('./00_rawdata')){
  dir.create('./00_rawdata')
}
setwd('./00_rawdata/')

### 常规情况（非高通量数据可以在GEO数据库直接用getGEO函数读取，否则要自己手动到网站下载表达矩阵）
library(GEOquery)
library(tidyverse)
library(lance)

GEO_data <- 'GSE119041'
gene_annotation <- 'GPL17692'

if(!dir.exists(paste0(GEO_data))){
  dir.create(paste0(GEO_data))
}
setwd(paste0(GEO_data))

gset <- getGEO(GEO_data , # 前面创建GEO_data对象
               destdir = '.',
               GSEMatrix = T,
               getGPL = F) # getGEO函数自动从官网中获取对应的数据集

expr <- exprs(gset[[1]]) # 将获取的数据集的表达矩阵提取出来，如果没有表达矩阵，说明是高通量数据需要到网站自己下载
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  expr[which(expr <= 0)] <- 0
  expr <- log2(expr + 1) 
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
expr <- as.data.frame(expr)
library(httr)
set_config(config(timeout = 300))
options(timeout = 300)

gpl <- getGEO(gene_annotation, destdir = '.')
gpl <- Table(gpl)
colnames(gpl)


gpl$gene_assignment<-data.frame(sapply(gpl$gene_assignment,function(x)unlist(strsplit(x," // "))[2]),
                                stringsAsFactors=F)[,1]
probe2symbol <- dplyr::select(gpl, 'ID', 'gene_assignment')
probe2symbol <- filter(probe2symbol, 'gene_assignment' != '')
probe2symbol <- separate(probe2symbol, 'gene_assignment', into = c('symbol', 'drop'), sep = '//')
probe2symbol <- dplyr::select(probe2symbol, -drop)
names(probe2symbol) <- c('ID', 'symbol')
probe2symbol <- probe2symbol[probe2symbol$symbol != ' --- ', ]
probe2symbol$symbol <- gsub(' ', '', probe2symbol$symbol)


dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)

dat <- dat %>%
  merge(probe2symbol, by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol, everything())%>%     ## 重新排列
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol, .keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

# protein_gene <- read.delim2('/data/nas2/database/gencode/PCG.xls(v36)')
# dat <- dat[rownames(dat) %in% protein_gene$gene_name, ]

a <- gset[[1]]
pd <- pData(a)
table(pd$title)
group <- data.frame(sample = pd$geo_accession, group = pd$title)
selected_columns <- group[grepl("^ESS",group$group),]
group <- selected_columns
table(group$group)
group$group <-ifelse(grepl("^Normal",group$group),"normal","Tumor")
table(group$group) 

dat <- dat[, colnames(dat) %in% group$sample]
dat <- na.omit(dat)
group <- group[group$sample %in% colnames(dat), ]

write.csv(dat, file = paste0('dat.', GEO_data, '.csv'))
write.csv(group, file = paste0('group.', GEO_data, '.csv'), row.names = F)

survival <- data.frame(sample = pd$geo_accession, 
                       OS = pd$characteristics_ch1.9,
                       OS.time= pd$characteristics_ch1.8)


survival <- as.data.frame(survival)
survival$OS <-ifelse(grepl("Alive$",survival$OS),"0","1")
survival$OS.time <- sub(".*:\\s*", "", survival$OS.time)
survival$OS.time <- as.numeric(survival$OS.time)
survival$OS.time <- survival$OS.time*30
survival <- na.omit(survival)
table(survival$OS)
# 0  1 
# 11 39
survival$OS<-as.numeric(survival$OS)
write.csv(survival, file = paste0('survival.', GEO_data, '.csv'), row.names = F)







