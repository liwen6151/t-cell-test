# 1.训练集------------------------
rm(list = ls())
setwd("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/")
if (! dir.exists("./04_risk_cox")){
  dir.create("./04_risk_cox")
}
setwd("./04_risk_cox")

library(GEOquery)
library(Biobase)
library(tidyverse)
library(lance)

group <- read.csv('../00_rawdata/TCGA/01.group_TCGA.UCEC.csv')
table(group$Type)
tumor_group <- group[group$Type == 'Tumor', ]

survival_data <- read.csv('../00_rawdata/TCGA/01.survival_TCGA.UCEC.csv', row.names = 1)
hub_gene <- read.csv('../02_hub_gene/01.gene_DGEs_TTKs_venn.csv', row.names = NULL)
dat <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F) %>% lc.tableToNum()
colnames(dat) <- gsub('.', '-', colnames(dat), fixed = T)
train_dat <- dat[rownames(dat) %in% hub_gene$symbol, ]
train_dat <- t(train_dat) # 转置
max(train_dat)
train_dat <- as.data.frame(train_dat)
train_dat$sample <- rownames(train_dat)
survival_data$sample <- rownames(survival_data)
train_dat <- merge(survival_data, train_dat, by='sample')
train_dat <- column_to_rownames(train_dat, var = 'sample')
train_dat <- train_dat[rownames(train_dat) %in% tumor_group$X, ]
train_data <- train_dat # 整理后的数据行名为样本名，第一列是OS状态，第二列是OS时间，从第三列开始就是取到的交集基因的表达


# 单因素COX回归分析-------------------------------------------------------------
## 单因素cox
library(survival)
library(survminer)

pfilter <- 0.01  # 设置p值的阈值
uniresult <- data.frame()   # 新建空白数据框
ph_res <- data.frame()
# i <- 2
# 使用for循环对输入数据中的基因依次进行单因素COX分析
# 单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(train_data[, 3:ncol(train_data)])){   
  unicox <- coxph(Surv(time = OS.time, event = OS) ~ train_data[, i], data = train_data)  
  cox_zph <- cox.zph(unicox)
  cox_table <- data.frame(cox_zph$table)
  if (cox_table$p[1] < 0.05) {
    next
  } else{
    rownames(cox_table) <- c(i, 'GLOBAL')
    ph_res <- rbind(ph_res, cox_table)
    unisum <- summary(unicox)   
    pvalue <- round(unisum$coefficients[, 5], 3) 
    if(pvalue < pfilter){ 
      uniresult <- rbind(uniresult,
                         cbind(gene = i,
                               HR = unisum$coefficients[, 2], # coef是公式中的回归系数b（有时候也叫beta值）, exp(coef)是cox模型中的风险比(HR), z代表wald统计量, 是coef除以其标准误se(coef)。
                               pvalue = unisum$coefficients[, 5], 
                               L95CI = unisum$conf.int[, 3], # lower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
                               H95CI = unisum$conf.int[, 4]
                         ))
    }
  }
  
}   
ph_res$symbol <- rownames(ph_res)
ph_res <- ph_res[ph_res$symbol %in% uniresult$gene, ]
uniresult$gene
#[1]"SPDEF"  "CENPF"  "CDKN2A" "E2F1"   "TRIP13" "AURKA"  "KLF2"   "TUBB2B" "CKMT1B" "SOX11"  "FSD1"   "ALK"    "OTC"   


uniresult <- column_to_rownames(uniresult, var = "gene")  ###看基因数目
uniresult <- lc.tableToNum(uniresult)
uniresult <- signif(uniresult, digits = 4)


res_results_0.05 <- as.data.frame(cbind(rownames(uniresult), uniresult$pvalue, paste0(uniresult$HR, " (", uniresult$L95CI, "-", uniresult$H95CI, ")"))) %>% lc.tableToNum()
names(res_results_0.05) <- c("gene", "p.value", "HR (95% CI for HR)")
res_results_0.05 <- column_to_rownames(res_results_0.05, var = "gene")

# 导出数据
write.csv(res_results_0.05, file = "univariate_cox_result.csv")
dim(res_results_0.05) # 4 2 行名是通过筛选的基因名（也就是所谓的具有预后价值的XX基因，第一列是其对应的P.value, 第二列是风险比及其置信区间）


res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)
res_results_0.05_2[, 1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]

hz <- paste(round(res_results_0.05_2$HR, 4),
            "(", round(res_results_0.05_2$HR.95L, 4),
            "-", round(res_results_0.05_2$HR.95H, 4),")", sep = "")


tabletext <- cbind(c(NA, "Gene", rownames(res_results_0.05_2)),
                   c(NA, "P value", ifelse(res_results_0.05_2$p.value<0.001,
                                           "< 0.001",
                                           round(res_results_0.05_2$p.value,3))),
                   c(NA, "Hazard Ratio(95% CI)", hz))

## 绘制森林图
library(forestplot)
p <- forestplot(labeltext = tabletext, 
                 mean=c(NA, NA, res_results_0.05_2$HR), 
                 lower=c(NA, NA, res_results_0.05_2$HR.95L), # 95%置信区间下限
                 upper=c(NA, NA, res_results_0.05_2$HR.95H), # 95%置信区间上限, 
                 col = fpColors(box="red", lines="royalblue", zero = "gray50"), 
                 graph.pos = 3,  # 为Pvalue箱线图所在的位置
                 is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
                 boxsize = 0.2, # 箱子大小
                 lwd.ci = 3,   # 线的宽度
                 ci.vertices.height = 0.08, #置信区间用线宽、高、型
                 ci.vertices = TRUE, 
                 zero = 1, # zero线宽 
                 lwd.zero = 0.5,    # 基准线的位置
                 colgap = unit(5, "mm"),    # 列间隙
                 xticks = c(0,1,2), # 横坐标刻度
                 lwd.xaxis = 2,  # X轴线宽
                 lineheight = unit(1.2,"cm"), # 固定行高
                 graphwidth = unit(.5,"npc"), # 图在表中的宽度比例
                 cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                 hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
                 # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                 #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                 #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
                 # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
                 #fpTxtGp函数中的cex参数设置各个组件的大小
                 txt_gp=fpTxtGp(label=gpar(cex=1),
                                ticks=gpar(cex=0.8, fontface = "bold"),
                                xlab=gpar(cex = 1, fontface = "bold"),
                                title=gpar(cex = 1.25, fontface = "bold")),
                 xlab="Hazard Ratio",
                 grid = T # 垂直于x轴的网格线，对应每个刻度
)
pdf(file = "01.univariate_cox_forest.pdf", height = 5, width = 9, onefile = F)
print(p)
dev.off()
png(filename = "01.univariate_cox_forest.png", height = 350, width = 600)
print(p)
dev.off()

ph_res
# chisq df         p symbol
# SPDEF  1.51690659  1 0.2180882  SPDEF
# CENPF  0.02691630  1 0.8696824  CENPF
# CDKN2A 0.13574816  1 0.7125450 CDKN2A
# E2F1   2.22682355  1 0.1356324   E2F1
# TRIP13 0.20963145  1 0.6470564 TRIP13
# AURKA  2.03231687  1 0.1539857  AURKA
# KLF2   0.01226168  1 0.9118285   KLF2
# TUBB2B 2.62041758  1 0.1054967 TUBB2B
# CKMT1B 1.72395605  1 0.1891846 CKMT1B
# SOX11  0.23786680  1 0.6257511  SOX11
# FSD1   0.50755806  1 0.4761979   FSD1
# ALK    0.01869602  1 0.8912415    ALK
# OTC    0.83814279  1 0.3599285    OTC

# 假定检验--------
save_plot<-function (p, filename, style, width = 5, height = 4) 
{
  if (style == "g") {
    ggsave(filename = paste0(filename, ".pdf"), p, width = width, 
           height = height)
    ggsave(filename = paste0(filename, ".png"), p, width = width, 
           height = height, dpi = 300)
  }
  else {
    pdf(file = paste0(filename, ".pdf"), width = width, family='Times',
        height = height)
    print(p, newpage = F)
    dev.off()
    png(filename = paste0(filename, ".png"), width = width, family='Times',
        height = height, units = "in", res = 300)
    print(p, newpage = F)
    dev.off()
  }
}

library(survminer)
library(GEOquery)
library(Biobase)
library(tidyverse)
library(ggsci)
library(data.table)
dir.create("./01.PH_plot")
res_results_0.05_2<- res_results_0.05_2 
gene_list <- rownames(res_results_0.05_2)
for(i in  rownames(res_results_0.05_2)){
  cox_formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", i))
  cox_model <- coxph(cox_formula, data = train_data)
  ph_test <- cox.zph(cox_model)
  print(ph_test)  # 查看检验结果
  cox_zph <- cox.zph(cox_model)
  p<-ggcoxzph(cox_zph)
  save_plot(p,filename = paste0("./01.PH_plot/",i),style = "x")
  ph_pvalues <- sapply(gene_list, function(gene) {
    cox_model <- coxph(as.formula(paste0("Surv(OS.time, OS) ~ ", gene)), data = train_data)
    cox.zph(cox_model)$table["GLOBAL", "p"]
  })
  # 保留 p ≥ 0.05（满足 PH 假设）
  res_results_ph <- res_results_0.05_2[ph_pvalues >= 0.05, ]
}


## 1.2 lasso模型----------------------
library(glmnet)
head(train_data, n=3) # 这个文件是前面整理好的的数据文件，行名为样本名，第一列是OS状态，第二列是OS时间，从第三列开始就是取到的交集基因名
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[, ph_res$symbol] # 选择前面筛选出的具有预后价值的基因所对应的样本
y_all <- subset(train_data, select = c(OS, OS.time))

## 拟合模型,glmnet只能接受矩阵形式的数据，因此要将数据框转成矩阵
fit <- glmnet(as.matrix(x_all), 
              Surv(y_all$OS.time, y_all$OS), 
              family = "cox")
print(fit)
fit$lambda


set.seed(300)
cvfit <- cv.glmnet(as.matrix(x_all),
                   Surv(y_all$OS.time, y_all$OS), 
                   nfolds=10,
                   family = "cox")

# 提取指定lambda时特征的系数
coef.min <- coef(cvfit, s = "lambda.min") ## lambda.min & lambda.1se 取一个

min = cvfit$lambda.min
cvfit$lambda.min # λ的最小值
# [1] 0.01254852
cvfit$lambda.1se # 最小值一个标准误的λ值
# [1] 0.08852731

# 找出那些回归系数没有被惩罚为0的
active.min <- which(coef.min@i != 0)
coef.min
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
# [1]  "SPDEF"  "CENPF"  "CDKN2A" "E2F1"   "KLF2"   "TUBB2B" "CKMT1B" "OTC" 


df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame

write.csv(df.coef, "Lasso_Coefficients.csv")
write.csv(lasso_geneids, './lasso_genes.csv')

### 绘图（优化·）
# 左边的虚线是均方误差最小时的λ，右边虚线是距离最小均方误差时的一个标准误的λ值
x <- coef(fit)
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1]
# extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
# compute L1 norm
tmp <- tmp[which(tmp$coef!="(Intercept)"),]
cvfit$lambda.min
head(tmp)
ggplot(tmp,aes(log(lambda),value,color = coef)) +
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype = 2)+
  geom_line(size=1) +
  xlab("Log (lambda)") +
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+
  # scale_color_manual(values= c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',
  #                              '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D',
  #                              ))+
  scale_x_continuous(expand = c(0.01,0.1))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=10,color='black'),
        legend.position = 'bottom')+
  annotate('text',x = -4,y= 0.85,
           label='Optimal Lambda = 0.01254852',
           color='black')+
  guides(col=guide_legend(ncol = 4))
ggsave("02.Lasso_model.pdf", height = 6, width = 5,family='Times')
ggsave("02.Lasso_model.png", height = 6, width = 5, units = 'in',dpi = 600) #字体不能改变，可将pdf文件转换


library(ggsci)
xx <- data.frame(lambda=cvfit[["lambda"]],cvm=cvfit[["cvm"]],cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],cvlo=cvfit[["cvlo"]],nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo,' vars')
ggplot(xx,aes(ll,cvm,color=NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=1)+
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log (lambda)")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+
  scale_color_manual(values = c(pal_npg()(10),pal_jco()(8)))+
  scale_x_continuous(expand = c(0.01,0.1))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'bottom')+
  annotate('text',x = -4,y=2.5,label = 'Optimal Lambda = 0.01254852' ,color='black')+
  guides(col=guide_legend(ncol = 4))
ggsave("03.Lasso_verify.pdf", height = 6, width = 5, family='Times')
ggsave("03.Lasso_verify.png", height = 6, width = 5, units = 'in', dpi = 600) #字体不能改变，可将pdf文件转换
