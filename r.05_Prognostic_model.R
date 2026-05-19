rm(list = ls())

setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./05_Prognostic_model")){
  dir.create("./05_Prognostic_model")
}
setwd("./05_Prognostic_model")


library(lance)
library(readxl)
library(readr)
library(tidyverse)


data <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F) %>% lc.tableToNum()
group <- read.csv('../00_rawdata/TCGA/01.group_TCGA.UCEC.csv')
table(group$Type)
tumor_group <- group[group$Type == 'Tumor', ]
survival_data <- read.csv('../00_rawdata/TCGA/01.survival_TCGA.UCEC.csv', row.names = 1)

train_dat <- data[,colnames(data) %in% tumor_group$X ]
data  <- train_dat
survival <- read.csv("../00_rawdata/TCGA/01.survival_TCGA.UCEC.csv")
colnames(survival)[colnames(survival) == "X"] <- "sample"
survival <- survival[survival$sample %in% tumor_group$X ,]
data <- data[,which(colnames(data)%in%survival$sample)]
gene<- read.csv("../04_risk_cox/lasso_genes.csv",header = T,row.names = 1)
rownames(gene) <- gene$x 
colnames(gene)
library(dplyr)  
## 合并生存数据
train_dat<-t(data[rownames(gene), colnames(data) %in% survival$sample]) %>% data.frame()
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_dat<-train_dat[,-1]
colnames(train_dat)


###验证集合------------
survival_test<-read.csv("../00_rawdata/GSE119041/survival.GSE119041.csv",header=TRUE, row.names=1)
survival_test$sample<-rownames(survival_test)
data_test <- read.csv('../00_rawdata/GSE119041/dat.GSE119041.csv',header=TRUE, row.names=1)
group <- read.csv('../00_rawdata/GSE119041/group.GSE119041.csv')
table(group$group)
survival_dat<-t(data_test)
survival_dat<-as.data.frame(survival_dat)
test_dat<-survival_dat[,rownames(gene)]

test_dat<-as.data.frame(test_dat)
test_dat$sample<-rownames(test_dat)
test_dat<-merge(survival_test,test_dat,by='sample')
rownames(test_dat)<-test_dat$sample
test_dat<-test_dat[,-1]
colnames(test_dat)

# 预处理数据
x_train <- subset(train_dat, select = -c(OS, OS.time))
colnames(x_train) <- gsub('\\.','-',colnames(x_train))
colnames(x_train)
y_train <- subset(train_dat, select = c(OS, OS.time))
colnames(y_train)

x_test <- subset(test_dat, select = -c(OS, OS.time))
colnames(x_test)
y_test <- subset(test_dat, select = c(OS, OS.time))
colnames(y_test)


library(randomForestSRC)
library(e1071)
library(survival)
library(caret)
library(survivalROC)
library(ggplot2)


# 设置随机种子
set.seed(1)

# 数据标准化
x_train <- scale(x_train)
x_test <- scale(x_test)

# 确保训练和测试数据的列顺序一致


# 将训练数据和测试数据合并
train_data <- data.frame(y_train, x_train)
test_data <- data.frame(y_test, x_test)
test_data<-na.omit(test_data)
table(test_data$OS)
save.image(file = 'train_data.RData')
load("train_data.RData")



# 重新训练模型                    
set.seed(1)
final_model <- rfsrc(Surv(OS.time, OS) ~ ., data = train_data,
                     ntree = 21, mtry =2)
final_model
# # # save.image(file = 'final_model.RData')
# # 绘制特征重要性图
# # 提取特征重要性
# importance <- vimp(final_model)$importance
# importance
# save.image(file = 'importance.RData')
# #
# # 转换为数据框
# importance_df <- data.frame(variable = names(importance), importance = importance)
# # 按重要性排序
# importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]
# 
# importance_plot <- ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
#   geom_bar(stat = 'identity', fill = 'blue') +
#   coord_flip() +
#   xlab('Features') +
#   ylab('Importance') +
#   ggtitle('Feature Importance Ranking')
# # 保存图像为 PDF
# ggsave("feature_importance1.pdf", plot = importance_plot, width = 8, height = 6)
# # 保存图像为 PNG
# ggsave("feature_importance1.png", plot = importance_plot, width = 8, height = 6)


# 获取训练集和测试集的风险评分
train_predictions <- predict(final_model, newdata = train_data)$predicted
test_predictions <- predict(final_model, newdata = test_data)$predicted


# 添加风险评分到训练集和测试集数据框中
train_dat$riskScore <- train_predictions
test_dat$riskScore <- test_predictions
identical(rownames(train_data), rownames(train_dat))

# 添加风险评分到训练集和测试集数据框中
train_data$riskScore <- train_predictions
test_data$riskScore <- test_predictions
write.csv(train_data,file = '01.train_risk.csv',row.names = T)
write.csv(test_data,file = '01.test_risk.csv',row.names = T)



###################test-----------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./05_Prognostic_model")){
  dir.create("./05_Prognostic_model")
}
setwd("./05_Prognostic_model")

if (! dir.exists("./test")){
  dir.create("./test")
}
setwd("./test")


Test_dat<-read.csv("../01.test_risk.csv",row.names = 1,check.names = F)

dat_score <- Test_dat[, c('OS', 'OS.time', 'riskScore')]
set.seed(1) 
library(survminer)
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS',variables = 'riskScore',minprop=0.1)
cutpoint <- res.cut$cutpoint[1]
cutpoint
dat_score$risk <- ifelse(dat_score$riskScore > res.cut$cutpoint$cutpoint, 1, 0)

# 18.49135
dat_score$risk <- ifelse(dat_score$risk == 1, 'high risk', 'low risk')
table(dat_score$risk)
write.csv(dat_score,file = 'risk.csv',row.names = T)
# KM曲线
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = dat_score)# 必要的库
library(survival)
library(survminer)


# 生成Kaplan-Meier曲线图和风险表
km_plot <- ggsurvplot(kmfit,
                      # pval = TRUE,
                      pval = surv_pvalue(kmfit,method = 'survdiff')[[4]],
                      conf.int = FALSE,
                      legend.labs = c("High risk", "Low risk"),
                      legend.title = "Risk score",
                      title = " ",
                      font.main = c(15, "bold"),
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      linetype = "strata",
                      ggtheme = theme_bw(),
                      palette = c("#f47e84","#6694e9"))
km_plot


# 使用ggpubr包中的ggarrange函数将两个图表合并
combined_plot <- ggarrange(km_plot$plot, km_plot$table, ncol = 1, heights = c(3, 1))

# 保存组合图表为PDF
ggsave("01.Test_KM.pdf", plot = combined_plot, width = 8, height = 6)

# 保存组合图表为PNG
ggsave("01.Test_KM.png", plot = combined_plot, width = 8, height = 6, dpi = 600)
library(timeROC)
library(survival)

# 构建timeroc
df<-dat_score
ROC <- timeROC(T=df$OS.time,
               delta=df$OS,
               marker=df$riskScore,
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c( 3*365, 5*365,7*365),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
ROC
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   #time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)

#添加标签信息
legend("bottomright",
       c(paste0("AUC at 3 years: ",round(ROC[["AUC"]][1],2)),
         paste0("AUC at 5 years: ",round(ROC[["AUC"]][2],2)),
         paste0("AUC at 7 years: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n")
pdf(file = "02.Test_ROC_Curve.pdf", width = 6, height = 6)
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   # time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    # add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 3 years: ", round(ROC[["AUC"]][1], 2)),
         paste0("AUC at 5 years: ", round(ROC[["AUC"]][2], 2)),
         paste0("AUC at 7 years: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n")
dev.off()

# 保存图形为 PNG 文件
png(file = "02.Test_ROC_Curve.png", width = 600, height = 600)
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   # time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    # add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 3 years: ", round(ROC[["AUC"]][1], 2)),
         paste0("AUC at 5 years: ", round(ROC[["AUC"]][2], 2)),
         paste0("AUC at 7 years: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n")
dev.off()

 
 
#风险评分分布--------
risk<-dat_score
id<-rownames(dat_score)
# 转换为因子，确保标签一致
risk$risk <- factor(risk$risk, levels = c("high risk", "low risk"), labels = c("High Risk", "Low Risk"))
library(ggplot2)
library(survminer)
riskScore <- risk$riskScore
risk$id <- rownames(risk)
risk_dis <- ggplot(risk, aes(x = reorder(id, riskScore),
                             y = riskScore,
                             color = factor(risk,
                                            levels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("High Risk" = "#f47e84", "Low Risk" = "#6694e9")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,2000,4000)],
                   labels = c(1,100,200,300,400,500,2000,4000),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  geom_hline(yintercept = cutpoint$cutpoint,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Risk Score",
       title = "Test") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  # 使用 legend.position = "bottom"
        legend.background = element_rect(color = "black", linewidth = 0.3),  # 使用 linewidth 替代 size
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(risk_dis)

# 生存状态分布
surv_stat <- ggplot(risk, aes(x = reorder(id, riskScore),
                              y = OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point(na.rm = F) +  # 忽略 NA 值
  scale_color_manual(values = c("#6694e9","#f47e84")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,2000,4000)],
                   labels = c(1,100,200,300,400,500,2000,4000),
                   expand = c(0.02, 0)) +
  ylim(c(0, max(risk$OS.time / 365, na.rm = TRUE))) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Survival time (Years)",
       title = "") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  # 使用 legend.position = "bottom"
        legend.background = element_rect(color = "black", linewidth = 0.3),  # 使用 linewidth 替代 size
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(surv_stat)


library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)
p
ggsave(filename = '03.Test_risk_and_survival_Distribution.pdf',p,w=6,h=7)
ggsave(filename = '03.Test_risk_and_survival_Distribution.png',p,w=6,h=7, units='in', dpi = 600)


risk2<-dat_score
group_risk <- risk2 %>% as.data.frame()
group_risk <- group_risk[order(group_risk$riskScore),]
colnames(group_risk)
group_risk2<-subset(group_risk, select = c(risk))
colnames(group_risk2) <- 'group'


#表达文件
dat <- Test_dat
gene<-read.csv("../../04_risk_cox/lasso_genes.csv")
coxGene<-gene$x
gene.expr <- dat[rownames(group_risk2),coxGene]
gene.expr<-as.data.frame(t(gene.expr))
mat <- gene.expr
mat <- t(scale(t(gene.expr)))  # 归一化,让数据分布更加均匀
# 给极大值和极小值赋值，让图片更好看
mat[mat < (-2)] <- -2
mat[mat > 2] <- 2


# 热图
library(ComplexHeatmap)
library(circlize)
pdf(file = '04.Test_heatmap.pdf',w=8,h=4,family='Times')
col_fun<-colorRamp2(c(-2,0,2),c("#7974a5", "white","#f4c466"))
HeatmapAnnotation(Group = group_risk2$group,
                  col = list(Group = c("high risk" = "#f4c466", "low risk" = "#7974a5")),
                  show_annotation_name = T,
                  show_legend = T) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),
          show_column_names = F,name = "expression",
          height = unit(6, "cm"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_names_side = "left",
          col = col_fun)
dev.off()

png(file = '04.Test_heatmap.png',w=8,h=4, units = 'in',res = 600,family='Times')
col_fun<-colorRamp2(c(-2,0,2),c("#7974a5", "white","#f4c466"))
HeatmapAnnotation(Group = group_risk2$group,
                  col = list(Group = c("high risk" = "#f4c466", "low risk" = "#7974a5")),
                  show_annotation_name = T,
                  show_legend = T) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),
          show_column_names = F,name = "expression",
          height = unit(6, "cm"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_names_side = "left",
          col = col_fun)
dev.off()





###################train-----------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/05_Prognostic_model/")
if (! dir.exists("./train")){
  dir.create("./train")
}
setwd("./train")


train_dat<-read.csv("../01.train_risk.csv",row.names = 1,check.names = F)

dat_score <- train_dat[, c('OS', 'OS.time', 'riskScore')]
set.seed(1) 
library(survminer)
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS',variables = 'riskScore',minprop = 0.2)
cutpoint <- res.cut$cutpoint[1]
cutpoint

# 19.00431
dat_score$risk <- ifelse(dat_score$riskScore > res.cut$cutpoint$cutpoint, 1, 0)
dat_score$risk <- ifelse(dat_score$risk == 1, 'high risk', 'low risk')
table(dat_score$risk)
write.csv(dat_score,file = 'risk.csv',row.names = T)

# KM曲线
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = dat_score)# 必要的库
library(survival)
library(survminer)

# 生成Kaplan-Meier曲线图和风险表
km_plot <- ggsurvplot(kmfit,
                      #pval = TRUE,
                      pval = surv_pvalue(kmfit,method = 'survdiff')[[4]],
                      conf.int = FALSE,
                      legend.labs = c("High risk", "Low risk"),
                      legend.title = "Risk score",
                      title = " ",
                      font.main = c(15, "bold"),
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      linetype = "strata",
                      ggtheme = theme_bw(),
                      palette = c("#f47e84","#6694e9"))
km_plot

# 使用ggpubr包中的ggarrange函数将两个图表合并
combined_plot <- ggarrange(km_plot$plot, km_plot$table, ncol = 1, heights = c(3, 1))

# 保存组合图表为PDF
ggsave("01.Train_KM.pdf", plot = combined_plot, width = 8, height = 6)

# 保存组合图表为PNG
ggsave("01.Train_KM.png", plot = combined_plot, width = 8, height = 6, dpi = 600)
library(timeROC)
library(survival)

# 构建timeroc
df<-dat_score
ROC <- timeROC(T=df$OS.time,
               delta=df$OS,
               marker=df$riskScore,
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c(3*365, 5*365, 7*365),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
ROC
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   #time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)

#添加标签信息
legend("bottomright",
       c(paste0("AUC at 3 years: ",round(ROC[["AUC"]][1],2)),
         paste0("AUC at 5 years: ",round(ROC[["AUC"]][2],2)),
         paste0("AUC at 7 years: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n")
pdf(file = "02.Train_ROC_Curve.pdf", width = 6, height = 6)
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   # time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    # add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 3 years: ", round(ROC[["AUC"]][1], 2)),
         paste0("AUC at 5 years: ", round(ROC[["AUC"]][2], 2)),
         paste0("AUC at 7 years: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n")
dev.off()

# 保存图形为 PNG 文件
png(file = "02.Train_ROC_Curve.png", width = 600, height = 600)
plot(ROC,
     time=3*365, col="red", lwd=2, title = "")   # time是时间点，col是线条颜色
plot(ROC,
     time=5*365, col="blue", add=TRUE, lwd=2)    # add指是否添加在上一张图中
plot(ROC,
     time=7*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 3 years: ", round(ROC[["AUC"]][1], 2)),
         paste0("AUC at 5 years: ", round(ROC[["AUC"]][2], 2)),
         paste0("AUC at 7 years: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n")
dev.off()



#风险评分分布--------
risk<-dat_score
id<-rownames(dat_score)
# 转换为因子，确保标签一致
risk$risk <- factor(risk$risk, levels = c("high risk", "low risk"), labels = c("High Risk", "Low Risk"))
library(ggplot2)
library(survminer)
riskScore <- risk$riskScore
risk$id <- rownames(risk)
risk_dis <- ggplot(risk, aes(x = reorder(id, riskScore),
                             y = riskScore,
                             color = factor(risk,
                                            levels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("High Risk" = "#f47e84", "Low Risk" = "#6694e9")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,2000,4000)],
                   labels = c(1,100,200,300,400,500,2000,4000),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  geom_hline(yintercept = cutpoint$cutpoint,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Risk Score",
       title = "Train") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  # 使用 legend.position = "bottom"
        legend.background = element_rect(color = "black", linewidth = 0.3),  # 使用 linewidth 替代 size
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(risk_dis)

# 生存状态分布
surv_stat <- ggplot(risk, aes(x = reorder(id, riskScore),
                              y = OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point(na.rm = F) +  # 忽略 NA 值
  scale_color_manual(values = c("#6694e9","#f47e84")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,2000,4000)],
                   labels = c(1,100,200,300,400,500,2000,4000),
                   expand = c(0.02, 0)) +
  ylim(c(0, max(risk$OS.time / 365, na.rm = TRUE))) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Survival time (Years)",
       title = "") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  # 使用 legend.position = "bottom"
        legend.background = element_rect(color = "black", linewidth = 0.3),  # 使用 linewidth 替代 size
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(surv_stat)


library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)
p

ggsave(filename = '03.Train_risk_and_survival_Distribution.pdf',p,w=6,h=7)
ggsave(filename = '03.Train_risk_and_survival_Distribution.png',p,w=6,h=7, units='in', dpi = 600)



# 热图
###########训练集预后基因表达-------------
risk2<-dat_score
group_risk <- risk2 %>% as.data.frame()
group_risk <- group_risk[order(group_risk$riskScore),]
colnames(group_risk)
group_risk2<-subset(group_risk, select = c(risk))
colnames(group_risk2) <- 'group'


#表达文件
dat <- train_dat
gene<-read.csv("../../04_risk_cox/lasso_genes.csv")
coxGene<-gene$x
gene.expr <- dat[rownames(group_risk2),coxGene]
gene.expr<-as.data.frame(t(gene.expr))
mat <- gene.expr
mat <- t(scale(t(gene.expr)))  # 归一化,让数据分布更加均匀
# 给极大值和极小值赋值，让图片更好看
mat[mat < (-2)] <- -2
mat[mat > 2] <- 2
library(ComplexHeatmap)
library(circlize)
pdf(file = '04.train_heatmap.pdf',w=8,h=4,family='Times')
col_fun<-colorRamp2(c(-2,0,2),c("#7974a5", "white","#f4c466"))
HeatmapAnnotation(Group = group_risk2$group,
                  col = list(Group = c("high risk" = "#f4c466", "low risk" = "#7974a5")),
                  show_annotation_name = T,
                  show_legend = T) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),
          show_column_names = F,name = "expression",
          height = unit(6, "cm"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_names_side = "left",
          col = col_fun)
dev.off()

png(file = '04.train_heatmap.png',w=8,h=4, units = 'in',res = 600,family='Times')
col_fun<-colorRamp2(c(-2,0,2),c("#7974a5", "white","#f4c466"))
HeatmapAnnotation(Group = group_risk2$group,
                  col = list(Group = c("high risk" = "#f4c466", "low risk" = "#7974a5")),
                  show_annotation_name = T,
                  show_legend = T) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),
          show_column_names = F,name = "expression",
          height = unit(6, "cm"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_names_side = "left",
          col = col_fun)
dev.off()
