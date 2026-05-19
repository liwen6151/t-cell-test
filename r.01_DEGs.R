#TCGA-UCEC-------------
rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if (! dir.exists("./01_DEGs/")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")

TCGA_data <- 'TCGA-UCEC'
threshold_p <- 0.05
threshold_logfc <- 2


if(!dir.exists(paste0(TCGA_data))){
  dir.create(paste0(TCGA_data))
}
setwd(paste0(TCGA_data))

library(DESeq2)

dat <- read.csv("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/00_rawdata/TCGA/01.count_TCGA.UCEC_mRNA.csv", row.names = 1)
colnames(dat) <- gsub('.', '-', colnames(dat), fixed = T)
colData <- read.csv("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/00_rawdata/TCGA/01.group_TCGA.UCEC.csv")
colnames(colData) <- c('sample', 'group')
table(colData$group)
# Normal  Tumor 
# 35       539 
colData <- colData[order(colData$group),]
dat <- dat[,colData$sample]
dim(dat) # [1] 19938   574

colnames(dat) == colData$sample
colData$group<-factor(colData$group)

rownames(colData) <- NULL
dds <- DESeqDataSetFromMatrix(countData=dat, 
                              colData=colData, 
                              design=~group)
dds_norm <- vst(dds)

# # # 提取标准化后的数据------
Controlized_counts <- assay(dds_norm)

dds$group <- relevel(dds$group, ref ="Normal")   #指定对照组
dds$group
dds <- DESeq(dds)
save.image("dds.Rdata")
load('./dds.Rdata')

# 提取差异结果------
res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
table(res$pvalue<0.05)
# FALSE  TRUE 
# 6192 13300  

DEG <- as.data.frame(res)
DEG <- na.omit(DEG)
dim(DEG)
head(DEG)
DEG <- cbind(symbol = rownames(DEG), DEG)
## 添加change列

logFC_cutoff <- 2

DEG$change = as.factor(
  ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
# DOWN   NOT    UP 
# 975 16901  1616

sig_diff <- subset(DEG,
                   DEG$pvalue<0.05 & abs(DEG$log2FoldChange) >  logFC_cutoff)
table(sig_diff$change)
# DOWN  NOT   UP 
# 975    0   1616

write.table(DEG, file = "01.DEG_all.csv", quote = F, sep = ",",row.names = FALSE)
write.table(sig_diff, file = "02.DEG_sig.csv", quote = F, sep = ",",row.names = FALSE)


# 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep <- DEG[rownames(DEG)%in%
                 rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],5),
                                head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],5))),]
volcano_plot <- ggplot(data = DEG, 
                       aes(x = log2FoldChange,
                           y = -log10(pvalue), 
                           color =change)) +
  scale_color_manual(values = c("#5F9EA0", "darkgray","#FF7F24")) +
  scale_x_continuous(breaks = c(-10,-5,-1,0,1,5,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-2,2),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,family = "Times" )+
  labs(x = "log2 (Fold Change)",
       y = "-log10 (pvalue)",title = "Tumor vs Normal",
       subtitle = paste(sprintf('P.value:  %.2f;', 0.05),
                        sprintf("log2FC: %.2f;", 2),
                        sprintf('Up: %1.0f; Down: %1.0f;', table(sig_diff$change)[3], table(sig_diff$change)[1]),
                        sprintf('Total: %1.0f', dim(sig_diff)[1])))
volcano_plot
ggsave('03.volcano.png', volcano_plot,width = 8, height = 6)
ggsave('03.volcano.pdf', volcano_plot,width = 8, height = 6)


# 密度热图------
library(ComplexHeatmap)

rt <- Controlized_counts
rt <- rt[rownames(sig_diff),]
group <- colData
dat_rep<-dat_rep[order(dat_rep$log2FoldChange),]
heat<-rt[rownames(dat_rep),]
x<-heat
mat <- t(scale(t(x)))#归一化
mat[mat < (-2)] <- (-2)
mat[mat > (2)] <- (2)

pdf('04.heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("Tumor" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#87CEEB", "white","#DB7093"))(100))
dev.off()

png('04.heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("Tumor" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#87CEEB", "white","#DB7093"))(100))
dev.off()
