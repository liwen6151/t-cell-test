rm(list = ls())
setwd("/data/nas1/lijia/58_KYGW-61101-6-NKT-KM112")

if (!dir.exists("03_nomo/")) {dir.create("03_nomo/")}
setwd("03_nomo/")

# 加载包
library(survival)
library(survminer)
library(timeROC)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(rms)
library(survcomp)


A <- read.csv("01.train_risk.csv")

#重命名为 sample
names(A)[names(A) == "X"] <- "sample"

B <- read.csv("phenotype.csv",row.names = 1, check.names = F)

B$stage <- case_when(
  grepl("Stage I[A-Z]?$", B$stage) ~ "Stage I",    # IA, IB, IC → I
  grepl("Stage II[A-Z]?$", B$stage) ~ "Stage II",  # II, IIA, IIB → II
  grepl("Stage III[A-Z]?[0-9]?$", B$stage) ~ "Stage III",  # III, IIIA...IIIC2 → III
  grepl("Stage IV", B$stage) ~ "Stage IV",         # IVB → IV
  TRUE ~ NA_character_
)



final_data <- merge(A, B,by = c("sample", "OS", "OS.time"), all = FALSE)
str(final_data)


# 精确提取你需要的列
final_data <- final_data[, c("sample", "OS", "OS.time", "riskScore", "age", "grade", "stage", "tumor_invasion")]


# ========================= 单因素 Cox 回归 =========================
clinic <- final_data   # 备份原始数据

# 将分类变量转换为因子并指定参考组
clinic$grade <- factor(clinic$grade, levels = c("G1", "G2", "G3", "High Grade"))  # 若存在 G4 可加入
clinic$stage <- factor(clinic$stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

# 1. riskScore (连续)
res.risk <- coxph(Surv(OS.time, OS) ~ riskScore, data = clinic) %>% summary
res.risk <- c(res.risk$conf.int[1, c(1,3,4)], res.risk$coefficients[1,5])  # HR, lower, upper, p

# 2. age (连续)
res.age <- coxph(Surv(OS.time, OS) ~ age, data = clinic) %>% summary
res.age <- c(res.age$conf.int[1, c(1,3,4)], res.age$coefficients[1,5])

# 3. grade (分类，G1为参考)
res.grade <- coxph(Surv(OS.time, OS) ~ grade, data = clinic) %>% summary
# 提取 G2 和 G3 的结果
grade_levels <- rownames(res.grade$conf.int)
res.grade2 <- c(res.grade$conf.int[1, c(1,3,4)], res.grade$coefficients[1,5])  # G2 vs G1
res.grade3 <- c(res.grade$conf.int[2, c(1,3,4)], res.grade$coefficients[2,5])  # G3 vs G1
res.grade4 <- c(res.grade$conf.int[3, c(1,3,4)], res.grade$coefficients[3,5])  # G3 vs G1


# 4. stage (分类，Stage I为参考)
res.stage <- coxph(Surv(OS.time, OS) ~ stage, data = clinic) %>% summary
stage_levels <- rownames(res.stage$conf.int)
res.stage2 <- c(res.stage$conf.int[1, c(1,3,4)], res.stage$coefficients[1,5])  # Stage II vs I
res.stage3 <- c(res.stage$conf.int[2, c(1,3,4)], res.stage$coefficients[2,5])  # Stage III vs I
res.stage4 <- c(res.stage$conf.int[3, c(1,3,4)], res.stage$coefficients[3,5])  # Stage IV vs I

# 5. tumor_invasion (连续)
res.invasion <- coxph(Surv(OS.time, OS) ~ tumor_invasion, data = clinic) %>% summary
res.invasion <- c(res.invasion$conf.int[1, c(1,3,4)], res.invasion$coefficients[1,5])

# 构建结果数据框（包含参考行）
res <- rbind(res.risk, res.age,
             res.grade2, res.grade3,res.grade4,
             res.stage2, res.stage3, res.stage4,
             res.invasion) %>% as.data.frame()

# 添加指标名称（注意与森林图显示一致）
res$Indicator <- c("Risk score", "Age",
                   "Grade (G2 vs G1)", "Grade (G3 vs G1)", "Grade (High Grade vs G1)",
                   "Stage (Stage II vs I)", "Stage (Stage III vs I)", "Stage (Stage IV vs I)",
                   "Tumor invasion")
colnames(res)[1:4] <- c("hr", "low", "up", "pv")

# 插入参考行（用于森林图分组标题）
ref_grade <- data.frame(hr = 1, low = 1, up = 1, pv = NA,
                        Indicator = "Grade\n(G1 Reference)")
ref_stage <- data.frame(hr = 1, low = 1, up = 1, pv = NA,
                        Indicator = "Stage\n(Stage I Reference)")
# 按顺序插入：在 grade 结果前插入 Grade 参考行，在 stage 结果前插入 Stage 参考行
library(dplyr)


res <- bind_rows(
  res[1:2, ],          # 1. Risk score + Age
  ref_grade,           # 2. Grade 分组标题
  res[3:5, ],          # 3. Grade 三个结果
  ref_stage,           # 4. Stage 分组标题
  res[6:8, ],          # 5. Stage 三个结果
  res[9, ]             # 6. Tumor invasion
)

res$p <- ifelse(is.na(res$pv), NA, 
                ifelse(res$pv < 0.001, "p < 0.001", paste0("p = ", signif(res$pv, 2))))
res$Indicator <- factor(res$Indicator, levels = rev(res$Indicator))

# 准备森林图数据
res2 <- data.frame(p.value = res$pv,
                   HR = res$hr,
                   HR.95L = res$low,
                   HR.95H = res$up,
                   Indicator = res$Indicator,
                   row.names = res$Indicator)

# 生成表格文本
hz <- paste(round(res2$HR,3), "(", round(res2$HR.95L,3), "-", round(res2$HR.95H,3), ")", sep = "")
tabletext <- cbind(c(NA, rownames(res2)),
                   c("P value", ifelse(res2$p.value < 0.001, "< 0.001", round(res2$p.value, 4))),
                   c("Hazard Ratio(95% CI)", hz))
tabletext[tabletext == "1(1-1)"] <- NA

# 设置汇总行（is.summary）：参考行需要加粗
is_summary <- c(T,                      # 表头
                rep(F, 2),              # risk, age
                T,                      # Grade 参考
                rep(F, 3),              # G2, G3
                T,                      # Stage 参考
                rep(F, 3),              # Stage II, III, IV
                F)                      # tumor_invasion


library(forestplot)
# 绘制森林图
pdf("02.univariate_cox_prog_forest.pdf", family = "Times", height = 12, width = 20, onefile = F)
forestplot(labeltext = tabletext, 
           graph.pos = 4,
           is.summary = is_summary,
           col = fpColors(box = "red", lines = "darkblue", zero = "gray50"),
           mean = c(NA, res2$HR),
           lower = c(NA, res2$HR.95L),
           upper = c(NA, res2$HR.95H),
           boxsize = 0.1, lwd.ci = 3,
           ci.vertices.height = 0.08, ci.vertices = TRUE,
           zero = 1, lwd.zero = 0.5,
           colgap = unit(5, "mm"),
           lwd.xaxis = 2,
           lineheight = unit(2.0, "cm"),
           graphwidth = unit(0.6, "npc"),
           cex = 1.2, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp = fpTxtGp(label = gpar(cex = 1.2, fontfamily = "Times"),
                            ticks = gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                            xlab = gpar(cex = 1.5, fontface = "bold", fontfamily = "Times"),
                            title = gpar(cex = 1.8, fontface = "bold", fontfamily = "Times")),
           xlab = "Hazard Ratio",
           grid = TRUE,
           title = "Univariate",
           clip = c(0, 30))
dev.off()

# 可选输出 PNG
png("02.univariate_cox_prog_forest.png", family = "Times", height = 12, width = 20, units = "in", res = 600)
forestplot(labeltext = tabletext, 
           graph.pos = 4,
           is.summary = is_summary,
           col = fpColors(box = "red", lines = "darkblue", zero = "gray50"),
           mean = c(NA, res2$HR),
           lower = c(NA, res2$HR.95L),
           upper = c(NA, res2$HR.95H),
           boxsize = 0.1, lwd.ci = 3,
           ci.vertices.height = 0.08, ci.vertices = TRUE,
           zero = 1, lwd.zero = 0.5,
           colgap = unit(5, "mm"),
           lwd.xaxis = 2,
           lineheight = unit(2.0, "cm"),
           graphwidth = unit(0.6, "npc"),
           cex = 1.2, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp = fpTxtGp(label = gpar(cex = 1.2, fontfamily = "Times"),
                            ticks = gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                            xlab = gpar(cex = 1.5, fontface = "bold", fontfamily = "Times"),
                            title = gpar(cex = 1.8, fontface = "bold", fontfamily = "Times")),
           xlab = "Hazard Ratio",
           grid = TRUE,
           title = "Univariate",
           clip = c(0, 30))
dev.off()

# ========================= 单因素 PH 检验 =========================
clinical <- na.omit(clinic)   # 删除缺失值

# 依次检验每个变量的比例风险假设
cox_risk <- coxph(Surv(OS.time, OS) ~ riskScore, data = clinical)
cox_risk_zph <- cox.zph(cox_risk, transform = "identity")
cox_table_risk <- t(as.data.frame(cox_risk_zph$table[-nrow(cox_risk_zph$table),]))

cox_age <- coxph(Surv(OS.time, OS) ~ age, data = clinical)
cox_age_zph <- cox.zph(cox_age, transform = "identity")
cox_table_age <- t(as.data.frame(cox_age_zph$table[-nrow(cox_age_zph$table),]))

cox_grade <- coxph(Surv(OS.time, OS) ~ grade, data = clinical)
cox_grade_zph <- cox.zph(cox_grade, transform = "identity")
cox_table_grade <- t(as.data.frame(cox_grade_zph$table[-nrow(cox_grade_zph$table),]))

cox_stage <- coxph(Surv(OS.time, OS) ~ stage, data = clinical)
cox_stage_zph <- cox.zph(cox_stage, transform = "identity")
cox_table_stage <- t(as.data.frame(cox_stage_zph$table[-nrow(cox_stage_zph$table),]))

cox_invasion <- coxph(Surv(OS.time, OS) ~ tumor_invasion, data = clinical)
cox_invasion_zph <- cox.zph(cox_invasion, transform = "identity")
cox_table_invasion <- t(as.data.frame(cox_invasion_zph$table[-nrow(cox_invasion_zph$table),]))

# 合并 PH 检验结果
library(dplyr)
cox_table_dan <- bind_rows(
  as.data.frame(cox_table_risk) %>% mutate(Variable = "riskScore"),
  as.data.frame(cox_table_age) %>% mutate(Variable = "age"),
  as.data.frame(cox_table_grade) %>% mutate(Variable = "grade"),
  as.data.frame(cox_table_stage) %>% mutate(Variable = "stage"),
  as.data.frame(cox_table_invasion) %>% mutate(Variable = "tumor_invasion")
)
write.csv(cox_table_dan, file = "PH_test_单.csv", row.names = FALSE)

# ========================= 多因素 Cox 回归 =========================
# 选择候选变量（可根据单因素结果取舍，这里以全部连续变量+分类变量为例）
features <- c("riskScore", "age", "grade", "tumor_invasion")
cox_data <- as.formula(paste0("Surv(OS.time, OS) ~ ", paste(features, collapse = "+")))
cox_multi <- coxph(cox_data, data = clinical)

# 多因素 PH 检验
cox_multi_zph <- cox.zph(cox_multi)
cox_table_multi <- cox_multi_zph$table[-nrow(cox_multi_zph$table), ]
write.csv(cox_table_multi, file = "PH_test_多.csv")

# 保留满足 PH 假设的变量 (p > 0.05)
valid_vars <- rownames(cox_table_multi)[cox_table_multi[, 3] > 0.05]
if (length(valid_vars) == 0) {
  message("No variable satisfies PH assumption, using all variables anyway.")
  valid_vars <- features
}
cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(valid_vars, collapse = "+")))
cox_multi_final <- coxph(cox_formula, data = clinical)

# 提取多因素结果
mul_res <- summary(cox_multi_final)
mul_coef <- mul_res$coefficients
mul_conf <- mul_res$conf.int

# 构建结果数据框（支持分类变量的多行）
multi_res <- data.frame(
  p.value = mul_coef[, 5],
  HR = mul_conf[, 1],
  HR.95L = mul_conf[, 3],
  HR.95H = mul_conf[, 4],
  Indicator = rownames(mul_coef)
)

# 如果 stage 或 grade 被包含，需要添加参考行（用于森林图分组）
ref_flag_grade <- any(grepl("grade", multi_res$Indicator))
ref_flag_stage <- any(grepl("stage", multi_res$Indicator))

# 为保持与单因素类似的显示风格，手动插入参考行（此处简化，直接使用所有结果）
# 更优雅的方法：直接将 Indicator 重命名为可读形式
multi_res$Indicator <- gsub("grade", "Grade ", multi_res$Indicator)
multi_res$Indicator <- gsub("stage", "Stage ", multi_res$Indicator)
multi_res$Indicator <- gsub("riskScore", "Risk score", multi_res$Indicator)
multi_res$Indicator <- gsub("age", "Age", multi_res$Indicator)

write.csv(multi_res, file = "03.multivariate_cox_prog_result.csv", quote = FALSE, row.names = FALSE)

# 绘制多因素森林图
hz_mul <- paste(round(multi_res$HR,3), "(", round(multi_res$HR.95L,3), "-", round(multi_res$HR.95H,3), ")", sep = "")
tabletext_mul <- cbind(c(NA, multi_res$Indicator),
                       c("P value", ifelse(multi_res$p.value < 0.001, "< 0.001", round(multi_res$p.value, 4))),
                       c("Hazard Ratio(95% CI)", hz_mul))

# 汇总行设置（若有多分类变量可手动添加参考行，这里简单处理，不加参考行）
is_summary_mul <- c(T, rep(F, nrow(multi_res)))

pdf("04.multivariate_cox_prog_forest.pdf", family = "Times", height = 8, width = 20, onefile = F)
forestplot(labeltext = tabletext_mul, 
           graph.pos = 4,
           is.summary = is_summary_mul,
           col = fpColors(box = "red", lines = "darkblue", zero = "gray50"),
           mean = c(NA, multi_res$HR),
           lower = c(NA, multi_res$HR.95L),
           upper = c(NA, multi_res$HR.95H),
           boxsize = 0.1, lwd.ci = 3,
           ci.vertices.height = 0.08, ci.vertices = TRUE,
           zero = 1, lwd.zero = 0.5,
           colgap = unit(5, "mm"),
           lwd.xaxis = 2,
           lineheight = unit(2.0, "cm"),
           graphwidth = unit(0.6, "npc"),
           cex = 1, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp = fpTxtGp(label = gpar(cex = 1.1, fontfamily = "Times"),
                            ticks = gpar(cex = 0.9, fontface = "bold", fontfamily = "Times"),
                            xlab = gpar(cex = 1.3, fontface = "bold", fontfamily = "Times"),
                            title = gpar(cex = 1.5, fontface = "bold", fontfamily = "Times")),
           xlab = "Hazard Ratio",
           grid = TRUE,
           title = "Multivariate",
           clip = c(0, 30))
dev.off()

png("04.multivariate_cox_prog_forest.png", family = "Times", height = 8, width = 20, units = "in", res = 600)
forestplot(labeltext = tabletext_mul, 
           graph.pos = 4,
           is.summary = is_summary_mul,
           col = fpColors(box = "red", lines = "darkblue", zero = "gray50"),
           mean = c(NA, multi_res$HR),
           lower = c(NA, multi_res$HR.95L),
           upper = c(NA, multi_res$HR.95H),
           boxsize = 0.1, lwd.ci = 3,
           ci.vertices.height = 0.08, ci.vertices = TRUE,
           zero = 1, lwd.zero = 0.5,
           colgap = unit(5, "mm"),
           lwd.xaxis = 2,
           lineheight = unit(2.0, "cm"),
           graphwidth = unit(0.6, "npc"),
           cex = 1, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp = fpTxtGp(label = gpar(cex = 1.1, fontfamily = "Times"),
                            ticks = gpar(cex = 0.9, fontface = "bold", fontfamily = "Times"),
                            xlab = gpar(cex = 1.3, fontface = "bold", fontfamily = "Times"),
                            title = gpar(cex = 1.5, fontface = "bold", fontfamily = "Times")),
           xlab = "Hazard Ratio",
           grid = TRUE,
           title = "Multivariate",
           clip = c(0, 30))
dev.off()



#列线图============================================================================
#clinical$sample <- rownames(clinical)
# table(clinical$pathologic_T)
# table(clinical$tumor_Stage)
df_all2 <- clinical[,c(1:4,8)]

# Nomogram
library(rms)
ddist <- datadist(df_all2)
options(datadist = "ddist")
#age
cox_data2 <- as.formula(paste0('Surv(OS.time, OS)~', paste(c("riskScore","tumor_invasion"), collapse = "+")))
res.cox <- psm(cox_data2, data = df_all2, dist = "lognormal")
surv <- Survival(res.cox)
function(x) surv(365*3, x)
function(x) surv(365*5, x)
function(x) surv(365*7, x)

nom.cox <- nomogram(res.cox, 
                    fun = list(function(x) surv(365*3, x), function(x) surv(365*5, x), function(x) surv(365*7, x)),
                    funlabel = c("3-year OS Probability", "5-year OS Probability", "7-year OS Probability"),
                    # maxscale = 10,    #想需要注意修改一下
                    lp = F,   #去除linear predictor
                    fun.at = c(0.01,seq(0.1,0.9,by=0.1),0.99))

png(filename = "05.nomogram_line_points.png", height = 7, width = 11,units = "in",res = 600)
par(family = "Times")
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.6)
dev.off()
pdf(file = "05.nomogram_line_points.pdf", height = 7, width = 11)
par(family = "Times")
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.6)
dev.off()



#校准曲线========================================================================
library(rms)
#features <- c("riskscore","Age","tumor_Stage","pathologic_T","pathologic_N","pathologic_M")
features <- c("riskScore","tumor_invasion")   #此处写哪些无所谓，在下一步ph检验那里会进行筛选

cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(features, collapse = "+")))
cox_more <- coxph(cox_data, data = clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]      #PH假定检验

cox_formula <- as.formula(paste("Surv(OS.time, OS)~",
                                paste(rownames(cox_table)[cox_table[,3]>0.05],    #PH假定检验p>0.05
                                      collapse = "+")))

ddist <- datadist(clinical)
options(datadist = "ddist")
cox_data <- cox_formula
coxm_3 <- cph(cox_data, data = clinical, surv = T, x = T, y = T, time.inc = 3*365)
cal_3 <- calibrate(coxm_3, u = 3*365, cmethod = "KM", m = 130, B = 130)
coxm_5 <- cph(cox_data, data = clinical, surv = T, x = T, y = T, time.inc = 5*365)
cal_5 <- calibrate(coxm_5, u = 5*365, cmethod = "KM", m = 130, B = 130)
coxm_7 <- cph(cox_data, data = clinical, surv = T, x = T, y = T, time.inc = 7*365)
cal_7 <- calibrate(coxm_7, u = 7*365, cmethod = "KM", m = 130, B = 130)

png(filename = "06.nomogram_predicted.png", family = "Times", height = 6.5, width = 6.5, units = "in", res = 600)
pdf(file = "06.nomogram_predicted.pdf", family = "Times", height = 6.5, width = 6.5)
par(mar=c(5,4,2,3),cex=1.5,family="Times")
plot(cal_3,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-7 year OS',#便签
     ylab='Actual 3-7 year OS (proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-7 year OS',#便签
     ylab='Actual 3-7 year OS (proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_7,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-7 year OS',#便签
     ylab='Actual 3-7 year OS (proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围



#加上图例
legend("bottomright", legend=c("3-year", "5-year", "7-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()


############### DCA决策曲线--------------

library(caret)
library(ggDCA)

cph1 <- rms::cph(Surv(OS.time, OS)~ riskScore, data = as.data.frame(df_all2))
cph2 <- rms::cph(Surv(OS.time, OS)~ tumor_invasion, data = as.data.frame(df_all2))

cph4 <- rms::cph(Surv(OS.time, OS)~ riskScore+tumor_invasion, data = as.data.frame(df_all2))


d_train <- dca(cph4,cph1,cph2, model.names = c("Nomogram","riskScore","Tumor_invasion"),times = "median")
p1 <- ggplot(d_train)+theme(legend.position = "top")
ggsave(paste0("07.dca.png"), p1, width = 8, height = 6, dpi = 300, units = "in", bg = "white")
ggsave(paste0("07.dca.pdf"), p1, width = 8, height = 6, units = "in", bg = "white")




##roc--------

library(timeROC)
os_risk_clinical <- df_all2
res.mul <- coxph(cox_data2, data = os_risk_clinical)
os_risk_clinical$pred <- predict(res.mul, newdata = os_risk_clinical, type = "lp")#计算列线图得分
cli_dat <- os_risk_clinical
cli_dat$OS.time <- cli_dat$OS.time/365
colnames(cli_dat)[colnames(cli_dat)=="pred"] <- "nomogram"
y1 <- 3
y2 <- 5
y3 <- 7
ROC.nomo <- timeROC(T=cli_dat$OS.time,
                    delta=cli_dat$OS,
                    marker=cli_dat$nomogram,
                    cause=1,
                    #weighting="marginal",
                    times=c(y1,y2,y3),
                    iid=TRUE)
# ROC.risk <- timeROC(T=cli_dat$DSS.time,
#                     delta=cli_dat$DSS,
#                     marker=cli_dat$riskScore,
#                     cause=1,
#                     weighting="marginal",
#                     times=c(y1,y2,y3),
#                     iid=TRUE)

# cli_dat$riskScore <- as.numeric(substr(cli_dat$riskScore,2,2))
# cli_dat$T.stage <- as.numeric(substr(cli_dat$T.stage,2,2))
# ROC.pn <- timeROC(T=cli_dat$DSS.time,
#                   delta=cli_dat$DSS,
#                   marker=cli_dat$N.stage,
#                   cause=1,
#                   weighting="marginal",
#                   times=c(y1,y2,y3),
#                   iid=TRUE)
# ROC.age <- timeROC(T=cli_dat$DSS.time,
#                   delta=cli_dat$DSS,
#                   marker=cli_dat$age,
#                   cause=1,
#                   weighting="marginal",
#                   times=c(y1,y2,y3),
#                   iid=TRUE)

#绘图
col=c("#0099DD", "#FF4858", "#FF9933")
j <- 2
for (i in c(y1,y2,y3)){
  pdf(file = paste0("06.Nomogram_",i,"_year_ROC.pdf"),width = 5,height = 5,family = "Times")
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("06.Nomogram_",i,"_year_ROC.png"),width= 5, height= 5, units="in", res=600,family = "Times")
  dev.control("enable")
  plot(ROC.nomo, time = i, col=col[j-1], lwd=2,title = "")
  #title(main = paste0(i,"-year survival"))
  # plot(ROC.age, time = i, col="#66C2A5", lwd=2, add = T)
  #plot(ROC.pn, time = i, col="orange", lwd=2, add = T)
  #plot(ROC.risk, time = i, col='#658fcb', lwd=2,title = "")
  #plot(ROC.pt, time = i, col="#FFAA33", lwd=2, add = T)
  #plot(ROC.nomo, time = i, col="#f38687", lwd=2, add = T)
  legend("bottomright",
         paste0("nomogram: AUC of ",i," year = ",sprintf("%.2f",ROC.nomo[["AUC"]][j-1])),
         # c(paste0("riskScore: AUC = ",sprintf("%.2f",ROC.risk[["AUC"]][j-1])),
         #   # paste0("Age: AUC = ",sprintf("%.2f",ROC.age[["AUC"]][j-1])),
         #   #paste0("Pathologic_N: AUC = ",sprintf("%.2f",ROC.pn[["AUC"]][j-1])),
         #   paste0("Pathologic_T: AUC = ",sprintf("%.2f",ROC.pt[["AUC"]][j-1])),
         #   paste0("nomogram: AUC = ",sprintf("%.2f",ROC.nomo[["AUC"]][j-1]))),
         #col=c("#658fcb","#FFAA33","#f38687"),
         col=col[j-1],
         lty=1, lwd=3,bty = "n")
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
  j <- j+1
}

