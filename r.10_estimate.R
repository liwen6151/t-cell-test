rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if(!dir.exists('./10_estimate')){
  dir.create('./10_estimate')
}
setwd('./10_estimate/')

# 加载包
library(estimate)
library(tidyverse)
library(tidyr)

risk <- read.csv("../05_Prognostic_model/train/risk.csv", header = T, row.names = 1)
risk$sample <- rownames(risk)
risk$sample <- gsub('.', '-', risk$sample, fixed = TRUE)
group <- dplyr::select(risk, sample, risk)
colnames(group) <- c('sample', 'group')
group$group <- ifelse(group$group == "high risk", 'high', 'low')
table(group$group)


dat <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', check.names = F, row.names = 1) %>% lc.tableToNum()
train_data <- dat
colnames(train_data) <- gsub('.', '-', colnames(train_data), fixed = T)


train_data <- train_data[, colnames(train_data) %in% group$sample]
identical(colnames(train_data), group$sample)

library(IOBR)

estimate <- deconvo_tme(eset = train_data, method = "estimate")
colnames(estimate)
estimate <- column_to_rownames(estimate, var = 'ID')
colnames(estimate) <- gsub('_estimate', '', colnames(estimate))


##estimate <- log2(estimate + 1)

estimate$sample <- rownames(estimate)

write.csv(estimate, 'estimate_res.csv')

temp1 <- merge(estimate, group, by = "sample")
temp1 <- temp1 %>% 
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "type",
    values_to = "scores"
  )

colnames(temp1)

library(rstatix) # Assuming 'add_significance' is part of the ggpubr package

# 差异分析
wilcox_res <- temp1 %>% 
  group_by(type) %>%
  wilcox_test(scores ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
## 导出数据
write.csv(wilcox_res, file = 'Estimate_wilcox_res.csv')



table(temp1$type)
## 为绘制小提琴图提供数据并整理
violin.cibersort <- temp1
# violin.cibersort <- dat.cibersort2[dat.cibersort2$check_points %in% DE.cibersort$check_points, ]
class(violin.cibersort$group) # 检查分组是否为因子，如果不是要转换成因子
violin.cibersort$group <- as.factor(violin.cibersort$group)

plot_data <- violin.cibersort[violin.cibersort$type != 'TumorPurity', ]
table(plot_data$group)

p1 <- ggplot(plot_data,aes(x = type, y = scores, fill = group)) +
  geom_violin(position = position_dodge(0.9),alpha = 0.5,
              width = 1,trim = T,
              color = NA) +
  geom_boxplot(width = .2,show.legend = F,
               position = position_dodge(0.9),
               color = 'grey20',alpha = 0.5,
               outlier.color = 'grey50') +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45,color = 'black', hjust = 1),#,hjust = 1
        legend.position = 'top',
  ) + 
  scale_fill_manual(values = c('high'='orange','low'='turquoise'),
                    name = '') +
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.format",
                     label.y = 6000,size = 5) +
  labs(x="")+
  theme(axis.title.x = element_text(color="black", size=14, face="bold")
  )
ggsave(filename = '01.ESTIMATEScore.pdf',p1,w=8,h=8)
ggsave(filename = '01.ESTIMATEScore.png',p1,w=8,h=8,dpi = 600)

