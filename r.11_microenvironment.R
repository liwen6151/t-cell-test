rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./11_microenvironment")){dir.create("./11_microenvironment")}
setwd("./11_microenvironment")

library(magrittr)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


group <- read.csv('../05_Prognostic_model/train/risk.csv', check.names = F)
group <- subset(group, select = c(1,5))
group$risk <- factor(ifelse(group$risk == "high risk", "High", "Low"))
colnames(group) <- c('sample', 'group')
train_data <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F)
max(train_data)
train_data <- train_data[, colnames(train_data) %in% group$sample]

gene_set <- read.table("/data/nas1/zhangzhaolei/pipline/microenvironment/mmc3.txt", header = T, sep ="\t")
table(gene_set$cells.type)
gene_list <- split(as.matrix(gene_set)[,1], gene_set[,2])

# 默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。
# 当输入表达式值是整数计数count时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
scores <-  gsva(as.matrix(train_data), gene_list, method = "ssgsea", kcdf='Gaussian', abs.ranking=TRUE)
rownames(scores)

write.csv(scores,
          file = "score.csv",
          quote = F)

# 热图----
colnames(group) <- c('sample', 'group')
group<-group[order(group$group),]

scores<-read.csv('./score.csv',check.names = F, row.names = 1)
scores<-scores[,group$sample]
annotation_col<-as.data.frame(group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(scores)

library(pheatmap)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(Group = c('Low'="#45a9b8",'High'="#f76a56"))
p <- pheatmap(
  scores,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
p
png(filename = "01.heatmap.png", height = 5, width = 8, units = 'in',res = 600,family='Times')
p
dev.off()
pdf(file = "01.heatmap.pdf", height = 5,width = 8,family='Times')
p
dev.off()


# 差异分析和箱线图
scores <- as.data.frame(t(scores))
scores$sample <- rownames(scores)
violin_dat <- merge(scores, group, by = 'sample') 
violin_dat <- violin_dat %>% 
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "ImmuneCell",
    values_to = "Score"
  )
colnames(violin_dat)
library(rstatix)

wilcox_res <- violin_dat %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

DE.wilcox <- wilcox_res[which(wilcox_res$p < 0.05), ] # 筛选差异显著的肿瘤细胞
# 18个
write.csv(wilcox_res, file = 'wilcox_res.csv')
write.csv(DE.wilcox, file = 'DE_wilcox_res.csv')

colnames(violin_dat)
## 为绘制小提琴图提供数据并整理
plot_data <- violin_dat[violin_dat$ImmuneCell %in% DE.wilcox$ImmuneCell, ]
class(plot_data$group) # 检查分组是否为因子，如果不是要转换成因子

## 绘制小提琴图
p <- ggplot(plot_data, aes(x = ImmuneCell, y = Score, fill = group)) +
  # geom_violin(trim=F, color="black", aes(fill = group)) + #绘制小提琴图, “color”设置小提琴图的轮廓线的颜色(不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值), 则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom = "errorbar",
               width = 0.5,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(x = ImmuneCell, y = Score, fill = group),
               width = 0.5,
               position = position_dodge(0.9), 
               outlier.shape = NA, 
               outlier.colour = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values = c('gold', "#355783"), name = "Group")+
  labs(title = "ssGSEA", x = "", y = "Score", size = 20) +
  stat_compare_means(data = plot_data,
                     mapping = aes(group = group),
                     label = "p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "black", face = "bold", size = 16), 
        axis.text.y = element_text(hjust = 0.5, colour ="black", face="bold", size=16), 
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 18),
        legend.title = element_text(face = "bold", size = 18),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p
ggsave(filename = '01.immue_cell_bar.pdf',p,w=12,h=8)
ggsave(filename = '01.immue_cell_bar.png',p,w=12,h=8,dpi = 600)


# 2. cor-------------------------------
cibersort_res <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/11_microenvironment/score.csv', row.names = 1, check.names = F)
cibersort_res <- t(cibersort_res) %>% as.data.frame()
DE_cell <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/11_microenvironment/DE_wilcox_res.csv', row.names = 1, check.names = F)
cibersort_res <- cibersort_res[, colnames(cibersort_res) %in% DE_cell$ImmuneCell]
cibersort_res$sample <- rownames(cibersort_res)

cor_data <- cibersort_res
rownames(cor_data) <- NULL
cor_data <- column_to_rownames(cor_data, var = 'sample')

# 然后使用column_to_rownames设置新的行名
cor_res <- round(cor(cor_data, method = 'spearman'), 3)
p.mat <- ggcorrplot::cor_pmat(cor_data)

write.csv(cor_res, 'cor_res.csv')
write.csv(p.mat, 'p.mat.csv')

library(corrplot)
library(RColorBrewer)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggExtra)

col1 = colorRampPalette(colors =c("blue","gray","red"),space="Lab") # space参数选择使用RGB或者CIE Lab颜色空间
pdf(file = "02.cor_heatmap.pdf",width = 8,height = 8,family="Times")
corrplot(corr = cor_res, 
         p.mat = p.mat,
         method = "pie",
         type = "upper", 
         tl.pos = "lt", tl.cex = 0.8, tl.col = "black", tl.srt = 45, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 0.7, 
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 0.7,
         col = col1(20))
dev.off()

png(file = "02.cor_heatmap.png",width = 8,height = 8,family="Times",units = 'in', res = 600)
corrplot(corr = cor_res, 
         p.mat = p.mat,
         method = "pie",
         type = "upper", 
         tl.pos = "lt", tl.cex = 0.8, tl.col = "black", tl.srt = 45, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 0.7, 
         col = col1(20)) # 用*作为显著性标签，sig.level参数表示0.05用*表示。0.01用**。0.001用***。pch.cex = 0.8,用于设置显著性标签的字符大小。
corrplot(corr = cor_res, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 0.7,
         col = col1(20))
dev.off()


#基因与细胞
hub_gene <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/04_risk_cox/lasso_genes.csv')
hub_data <- train_data[rownames(train_data) %in% hub_gene$x, ] %>% t() %>% as.data.frame()

expr <-hub_data
tiics_result <- read.csv('./score.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr),]

diff <- read.csv('./DE_wilcox_res.csv')
tiics_result <- tiics_result[,diff$ImmuneCell]

identical(rownames(expr), rownames(tiics_result))

cor_r <- cor(expr,tiics_result,method = "spearman") 
cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(expr)))

library(psych)
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
write.csv(cor_dat,"03.correlation_cor.csv")

# 相关性热图带显著性----
data <- read.csv('03.correlation_cor.csv',row.names = 1,check.names = F)
data <- data %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    Pvalue <= 0.001 ~ "\n***", #P<0.001就显示回车加三个星号
    between(Pvalue, 0.001, 0.01) ~ "\n**", #P为0.001-0.01 显示回车加两个*号
    between(Pvalue, 0.01, 0.05) ~ "\n*",  #P为0.01-0.05 显示回车加一个星号
    T ~ ""))
p <-
  ggplot(data, aes(gene, cell)) +
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
  
  geom_text(aes(label = text),col ="black",size = 5) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(hjust = 1, size = 10, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 10, face = "bold") ,#调整y轴文字
        legend.title=element_text(size=15,family = "Times", face = "bold") , 
        text=element_text(family = 'Times'),
        legend.text=element_text(size=15,family = "Times", face = "bold")) + # 使用简洁的主题
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n"," *** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") #
p
ggsave(file=paste0('03.correlation_biomarker.png'), height = 8, width = 10, p)
ggsave(file=paste0('03.correlation_biomarker.pdf'), height = 8, width = 10, p)



