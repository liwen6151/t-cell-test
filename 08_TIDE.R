rm(list = ls())
setwd("/data/nas1/lijia/58_KYGW-61101-6-NKT-KM112")

if (!dir.exists("05_TIDE/")) {dir.create("05_TIDE/")}
setwd("05_TIDE/")

library(rstatix)
library(magrittr)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


group <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/05_Prognostic_model/train/risk.csv', check.names = F)
group <- subset(group, select = c(1,5))
group$risk <- factor(ifelse(group$risk == "high risk", "High", "Low"))
colnames(group) <- c('sample', 'group')
train_data <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F)
max(train_data)
train_data <- train_data[, colnames(train_data) %in% group$sample]


# 3. 计算CYT评分
TIDE_score <- read.csv('TIDE.csv')


res.cibersort <- TIDE_score


# 画图----
library(RColorBrewer)

res.cibersort2 <- res.cibersort


# 画图----
dat.cibersort <- res.cibersort2

dat.cibersort <- merge(group, dat.cibersort, by = "sample")


dat.cibersort2 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

colnames(dat.cibersort2)
stat_cibersort <- dat.cibersort2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.csv(stat_cibersort,file = '01.stat.cibersort.csv',row.names = T)
DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]  # 11 10
write.csv(DE.cibersort,file = '02.DE.cibersort.csv',row.names = T)

sig_cells <- DE.cibersort$ImmuneCell
# 1. 定义核心数据（re1：样本+分组+免疫细胞得分）
re1 <- dat.cibersort  # 复用整理后的数据
re1 <- re1[, c("sample", "group", sig_cells)]
# 2. 定义绘图数据（box：长格式，列名适配原逻辑）
box <- dat.cibersort2 %>%
  rename(
    Cell_type = ImmuneCell,  # 免疫细胞类型
    Proportion = Score       # 比例值
  ) %>%
  dplyr::filter(Cell_type %in% sig_cells) %>%
  mutate(Cell_type = factor(Cell_type, levels = colnames(re1)[-c(1:2)]))  # 按免疫细胞排序

# 3. 统计检验与颜色映射（以Low为参照组）
# 计算Low组样本数量（关键：替换原Low为Low）
Low_number <- table(re1$group)[["Low"]]
if (is.na(Low_number)) stop("未找到'Low'分组！请检查group$group的实际值")

d <- c()  # x轴标签颜色
e <- c()  # 箱线图边框颜色
# 遍历所有免疫细胞类型（排除sample和group列）
for (i in colnames(re1)[-c(1:2)]) {
  # 提取Low组和High组的得分（稳健提取，不依赖样本顺序）
  Low_scores <- re1[re1$group == "Low", i]  
  High_scores <- re1[re1$group == "High", i]  
  
  # 计算两组中位数
  median_Low <- median(Low_scores)
  median_High <- median(High_scores)
  
  # wilcox检验（Low vs High）
  t_res <- wilcox.test(Low_scores, High_scores)
  p_val <- t_res$p.value
  
  # 颜色规则：无差异→黑色；Low高（显著）→绿色；High高（显著）→紫色
  border_color <- if_else(
    p_val > 0.05, 
    "black", 
    if_else(median_Low > median_High, "#68A180", "#8C549C")
  )
  d <- c(d, border_color)                # x轴标签颜色（1个免疫细胞1种颜色）
  e <- c(e, border_color, border_color)  # 箱线图边框颜色（1个免疫细胞对应2组，需2种颜色）
}

# 4. 先定义ggplot对象（避免重复代码）
plot_immune_diff <- ggplot(box, aes(x = Cell_type, y = Proportion, fill = group)) +
  geom_boxplot(alpha = 0.7, col = e) +  # 箱线图（边框颜色反映差异方向）
  scale_fill_manual(values = c("#53A85F", "#E95C59"), name = "Group") +  # Low绿、High红
  labs(x = "", y = "CYT Score") +  # 坐标轴标签
  theme_bw() +
  theme(
    legend.position = 'top',          # 图例在顶部
    text = element_text(size = 20),   # 全局文本大小
    axis.line = element_line(color = "black"),  # 坐标轴线条
    axis.title = element_text(face = "bold", size = 22),  # 坐标轴标题（加粗）
    axis.text.x = element_text(       # x轴标签（免疫细胞类型）
      size = 16, vjust = 1, hjust = 1, angle = 0, 
      colour = d, face = "bold"       # 颜色与差异方向对应
    ),
    axis.text.y = element_text(size = 12),  # y轴标签大小
    panel.border = element_blank(),         # 移除面板边框
    panel.background = element_blank(),     # 移除面板背景
    panel.grid.major = element_blank(),     # 移除主网格线
    panel.grid.minor = element_blank()      # 移除次网格线
  ) +
  stat_compare_means(                    # 添加显著性标记
    label = "p.signif",                  # 显示*符号（ns/~/*/**/***）
    method = "wilcox.test",              # 与前面一致的检验方法
    hide.ns = TRUE,                      # 隐藏不显著（p>0.05）的标记
    size = 6                             # 符号大小
  )

# 5. 保存为PDF格式（矢量图，适合论文排版）
ggsave(
  filename = "03.cibersort.plot-2.pdf",
  plot = plot_immune_diff,
  width = 6,  # 图宽（与PNG一致，保证比例）
  height = 6,  # 图高
  device = "pdf"  # 明确指定设备为PDF
)

# 6. 保存为PNG格式（位图，适合快速预览、PPT展示）
ggsave(
  filename = "03.cibersort.plot-2.png",
  plot = plot_immune_diff,
  width = 6,  # 与PDF保持相同宽高比，避免拉伸
  height = 6,
  device = "png",  # 明确指定设备为PNG
  dpi = 300  # 分辨率（300dpi适合大多数场景，可根据需求调整）
)



# 免疫细胞相关性热图----
# 细胞相关性
tiics_result <- res.cibersort

rownames(tiics_result) <- tiics_result$sample
tiics_result$sample <- NULL   # 删除多余的列



library(psych)

hub_gene <- read.csv('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/04_risk_cox/lasso_genes.csv')
hub_gene <- train_data[rownames(train_data) %in% hub_gene$x, ] %>% t() %>% as.data.frame()
hub_gene$X <- rownames(hub_gene)

riskScore <- read.csv("/data/nas1/lijia/58_KYGW-61101-6-NKT-KM112/03_nomo/01.train_risk.csv")
riskScore <- riskScore[, c("X","riskScore")]


expr_Core<- merge(
  x = riskScore,        # 包含 X = 样本ID
  y = hub_gene,        # 包含 sampleID = 样本ID
  by = "X",           # 第一个数据的匹配列
  all = TRUE            # 保留所有样本
)

rownames(expr_Core) <- expr_Core$X
expr_Core$X <- NULL   # 删除多余的X列

# 可选：查看数据结构
str(tiics_result)
str(expr_Core)



# 3. 计算相关性
# 使用psych包的corr.test函数计算两组变量之间的相关性
# method = "pearson" 指定使用皮尔逊相关系数
# adjust = "fdr" 指定使用FDR方法对p值进行多重检验校正
data.corr <- corr.test(tiics_result, expr_Core, method = "spearman", adjust = "fdr")

# 提取相关系数矩阵和校正后的p值矩阵
correlation_matrix <- data.corr$r  # 相关系数
p_value_matrix <- data.corr$p      # 校正后的p值

# 4. 整理结果为长格式表格
# 将相关系数矩阵转换为长格式
cor_r_long <- correlation_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Immcell") %>%  # 将行名（Immcell基因）变为一列
  tidyr::gather(., key = "Core_Gene", value = "Correlation", -Immcell)  # 宽格式转长格式

# 将p值矩阵转换为长格式
cor_p_long <- p_value_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Immcell") %>%
  tidyr::gather(., key = "Core_Gene", value = "P_value", -Immcell)

# 合并相关系数和p值，并按指定顺序排列列
# 使用dplyr的inner_join更安全，确保只合并匹配的行
cor_dat <- inner_join(cor_r_long, cor_p_long, by = c("Immcell", "Core_Gene")) %>%
  select(Immcell, Core_Gene, Correlation, P_value) # 重新排列列的顺序

# 5. 保存结果到CSV文件
write.csv(cor_dat, file = "05.correlation_results.csv", row.names = FALSE)



rm(list = ls())


#library(corrplot)
library(RColorBrewer)
library(ggcorrplot)
library(ggplot2)
library(ggpubr)
library(ggExtra)


# 相关性热图带显著性----
data <- read.csv('05.correlation_results.csv',check.names = F)
# 新增：保留两位小数的相关系数，并与星号拼接
data <- data %>%
  mutate(
    corr_label = sprintf("%.2f", Correlation),  # 格式化相关系数为两位小数
    text = case_when(  # 设置显著性星号
      P_value <= 0.001 ~ "***",
      between(P_value, 0.001, 0.01) ~ "**",
      between(P_value, 0.01, 0.05) ~ "*",
      TRUE ~ ""
    ),
    # 关键修改：仅对有显著性的结果显示标签，无显著性则为空
    full_label = ifelse(text != "", paste0(corr_label, "\n", text), "")
  )
p <- ggplot(data, aes(x = Immcell, y = Core_Gene)) +
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  
  # 关键修改：显示合并后的标签（相关系数+星号）
  geom_text(aes(label = full_label), col = "black", size = 4) +  # 可根据需要调整字体大小
  
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 15, family = "Times", face = "bold"),
    text = element_text(family = 'Times'),
    legend.text = element_text(size = 15, family = "Times", face = "bold")
  ) +
  labs(fill = paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n"," *** p < 0.001","\n\n","Correlation")) +
  scale_x_discrete(position = "top")

p

ggsave(file=paste0('06.correlation_biomarker.png'), height = 8, width =4, p)
ggsave(file=paste0('06.correlation_biomarker.pdf'), height =8, width = 4, p)
