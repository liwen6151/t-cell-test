#  相关性图--------------------------------------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./14_cor/")){
  dir.create("./14_cor")
}
setwd("./14_cor")

library(tidyverse)
library(lance)
library(psych)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(circlize)
library(ggcorrplot)


train_data <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F)
max(train_data)

hub_gene <- read.csv('../04_risk_cox/lasso_genes.csv')
hub_gene <- hub_gene$x
hub_data <- train_data[rownames(train_data) %in% hub_gene, ] %>% t() %>% as.data.frame()
hub_data$sample <- rownames(hub_data)

cor_data <- hub_data
cor_data <- cor_data[,-9]
cor_data <- as.matrix(cor_data)

cor_res <- round(cor(cor_data, method = 'spearman'), 3)
p.mat <- ggcorrplot::cor_pmat(cor_data)
write.csv(p.mat,file="01.correlation_pvalue.csv",quote=F)
write.csv(cor_res, './02.correlation_relation.csv', row.names = F)

## 相关性热图——获取下三角形
library(corrplot)
library(RColorBrewer)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggExtra)
col1 = colorRampPalette(colors =c("blue","gray","red"),space="Lab") # space参数选择使用RGB或者CIE Lab颜色空间

pdf(file = "03.cor_heatmap.pdf",width = 7,height = 7,family="Times")
corrplot(corr = cor_res,
         p.mat = p.mat,
         method = "pie",
         type = "upper",
         tl.pos = "lt", tl.cex = 1, tl.col = "black", tl.srt = 45, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.2,
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res,
         type="lower",
         add=TRUE,
         method="number",
         tl.pos = "n",
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1.2,
         col = col1(20))
dev.off()

png(file = "03.cor_heatmap.png",width = 7,height =7,family="Times",units = 'in', res = 600)
corrplot(corr = cor_res,
         p.mat = p.mat,
         method = "pie",
         type = "upper",
         tl.pos = "lt", tl.cex = 1, tl.col = "black", tl.srt = 45, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.2,
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res,
         type="lower",
         add=TRUE,
         method="number",
         tl.pos = "n",
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1.2,
         col = col1(20))
dev.off()


#functional similarity----------------------------------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./14_cor/")){
  dir.create("./14_cor")
}
setwd("./14_cor")

library(GOSemSim)
library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(ggplot2)
library(org.Hs.eg.db)

GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库

inter_gene <- read.csv('../04_risk_cox/lasso_genes.csv', row.names = 1)

## gene ID转换
gene <- bitr(inter_gene$x, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)


bp <- godata(OrgDb = GO_database, ont = 'BP', computeIC = FALSE)
cc <- godata(OrgDb = GO_database, ont = 'CC', computeIC = FALSE)
mf <- godata(OrgDb = GO_database, ont = 'MF', computeIC = FALSE)

simbp <- mgeneSim(gene$ENTREZID,
                  semData = bp,
                  measure = 'Wang',
                  drop = NULL, 
                  combine = 'BMA')
simcc <- mgeneSim(gene$ENTREZID,
                  semData = cc,
                  measure = 'Wang',
                  drop = NULL, 
                  combine = 'BMA')
simmf <- mgeneSim(gene$ENTREZID,
                  semData = mf,
                  measure = 'Wang',
                  drop = NULL, 
                  combine = 'BMA')
fsim <- (simbp * simcc * simmf)^(1/3)   #simbp *

colnames(fsim) <- gene$SYMBOL
rownames(fsim) <- gene$SYMBOL

for (i in 1 : ncol(fsim)) {
  fsim[i,i] <- NA
}

dat <- melt(fsim)
dat <- dat[!is.na(dat$value), ]
write.csv(dat, './04.functional_correlation_relation.csv', row.names = F)
dat <- dat[, c(1, 3)]

dat_mean <- aggregate(value ~ Var1, dat , mean)
m <- dat_mean$value
names(m) <- dat_mean$Var1

dat$Var1 <- as.factor(dat$Var1)

p1 <- ggplot(dat, aes(x = Var1, y = value, fill = factor(Var1)))+
  scale_fill_brewer(palette = 'Set2')+
  geom_boxplot()+
  coord_flip()+
  labs(x = '', y = '', title = 'Functional similarity between genes')+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face="bold",color="black",family = "Times",size=12),
    plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
    axis.text.x = element_text(face = "bold",color = "black",size = 12),
    axis.text.y = element_text(face = "bold",color = "black",size = 12),
    axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
    axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
    plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 12, face = "italic", colour = "black"),
  )
p1
ggsave(filename = '04.GOSemSim.pdf', p1, w=8, h=8)
ggsave(filename = '04.GOSemSim.png', p1, w=8, h=8, dpi = 600)
