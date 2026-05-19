# 1. 药物IC50----------------------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./09_Drug")){
  dir.create("./09_Drug")
}
setwd("./09_Drug")

library(tidyverse)
library(lance)
library(oncoPredict)
library(ggplot2)

GDSC2_expr <- readRDS('/data/nas1/zhangzhaolei/pipline/IC50/DataFiles/Training_Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_res <- readRDS('/data/nas1/zhangzhaolei/pipline/IC50/DataFiles/Training_Data/GDSC2_Res.rds')
dat <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv',check.names = F, row.names = 1)

train_data <- dat


colnames(train_data) <- gsub('.', '-', colnames(train_data), fixed = TRUE)
risk <- read.csv('../05_Prognostic_model/train/risk.csv',  check.names = F) #%>% dplyr::select(id, riskScore)
risk <- risk[,c(1,4)]
colnames(risk) <- c('sample', 'riskscore')
train_data <- train_data[, colnames(train_data) %in% risk$sample]

calcPhenotype(trainingExprData = GDSC2_expr, 
              trainingPtype = GDSC2_res, 
              testExprData = as.matrix(train_data), 
              batchCorrect = 'eb', 
              powerTransformPhenotype = TRUE, 
              removeLowVaryingGenes = 0.2, 
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' 
)


## 1.1 药物IC50高低风险组wilcox检验-----------------------------------------------------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./09_Drug")){
  dir.create("./09_Drug")
}
setwd("./09_Drug")

library(tidyverse)
library(lance)
library(oncoPredict)
library(ggplot2)


IC50_res <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', row.names = 1, check.names = FALSE)
IC50_res$sample <- rownames(IC50_res)

risk <- read.csv('../05_Prognostic_model/train/risk.csv',  check.names = F) #%>% dplyr::select(id, riskScore)
risk <- risk[,c(1,5)]
colnames(risk) <- c('sample', 'group')
# risk$group <- ifelse(risk$group == 0, 'high', 'low')
table(risk$group)
risk$group <- as.factor(risk$group)

dat.IC50 <- merge(IC50_res, risk, by = "sample")
dat.IC50_2 <- dat.IC50 %>%
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "drug",
    values_to = "Score"
  )

colnames(dat.IC50_2)

library(rstatix)
## 差异分析
stat_res <- dat.IC50_2 %>%
  group_by(drug) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

DE.res <- stat_res[which(stat_res$p < 0.05), ]

## 导出数据
write.csv(stat_res, file = 'stat.IC50.csv')
write.csv(DE.res, file = 'DE.IC50.csv')

colnames(DE.res)

DE.res <- DE.res[order(DE.res$p), ]
DE.res <- DE.res[1 : 20, ]

## DE_1
## 为绘制小提琴图提供数据并整理
# violin.cibersort1 <- dat.IC50_2
violin.cibersort1 <- dat.IC50_2[dat.IC50_2$drug %in% DE.res$drug, ]
violin.cibersort1 <- separate(violin.cibersort1, drug, into = c('drug', 'id'), sep = '_')
class(violin.cibersort1$group) # 检查分组是否为因子，如果不是要转换成因子

library(ggpubr)
## 绘制小提琴图
p1 <- ggplot(violin.cibersort1, aes(x = drug, y = Score, fill = group)) +
  # geom_violin(trim=F, color="black", aes(fill = group)) + #绘制小提琴图, “color”设置小提琴图的轮廓线的颜色(不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值), 则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom = "errorbar",
               width = 0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(x = drug, y = Score, fill = group),
               width = 0.2,
               position = position_dodge(0.9),
               outlier.shape = NA,
               outlier.colour = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values = c('gold', "#355783"), name = "Group")+
  labs(title = "IC50", x = "", y = "IC50", size = 20) +
  stat_compare_means(data = violin.cibersort1,
                     mapping = aes(group = group),
                     label = "p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "black", face = "bold", size = 10),
        axis.text.y = element_text(hjust = 0.5, colour ="black", face="bold", size=12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1
ggsave(filename = '01.IC50.pdf',p1,w=12,h=8)
ggsave(filename = '01.IC50.png',p1,w=12,h=8,dpi = 600)

