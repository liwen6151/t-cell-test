rm(list = ls())
setwd("/data/nas1/lijia/58_KYGW-61101-6-NKT-KM112")

if (!dir.exists("06_正态性检验和方差齐性检验/")) {dir.create("06_正态性检验和方差齐性检验/")}
setwd("06_正态性检验和方差齐性检验/")


# 1. 录入数据
# SPDEF
SPDEF_C <- c(0.942682752,0.954691221,1.132991853,1.050604248,0.919029926)
SPDEF_U <- c(1.182943271,1.227403386,1.292882879,1.07254862,1.335938059)

# CENPF
CENPF_C <- c(0.981486052,0.994218646,1.059037371,0.995606642,0.969651288)
CENPF_U <- c(1.000695423,1.057949833,1.30827587,1.107519324,1.124727116)

# CDKN2A
CDKN2A_C <- c(1.082723382,0.885031026,0.881823105,0.882340226,1.268082261)
CDKN2A_U <- c(1.177983854,1.003619284,0.936506515,1.039271396,1.363618962)

# E2F1
E2F1_C <- c(0.947360448,0.891899868,1.022410133,0.942612476,1.195717075)
E2F1_U <- c(1.070144674,1.260329761,1.083063529,1.211922882,1.297974981)

# KLF2
KLF2_C <- c(1.0830909,0.963966005,1.074347987,0.968165132,0.910429977)
KLF2_U <- c(0.692226314,0.607282024,0.770283685,0.68205064,0.779876053)

# 2. 正态性检验 Shapiro-Wilk
shapiro.test(SPDEF_C); shapiro.test(SPDEF_U)
shapiro.test(CENPF_C); shapiro.test(CENPF_U)
shapiro.test(CDKN2A_C); shapiro.test(CDKN2A_U)
shapiro.test(E2F1_C); shapiro.test(E2F1_U)
shapiro.test(KLF2_C); shapiro.test(KLF2_U)

# 3. 方差齐性检验 Levene
library(car)
leveneTest(c(SPDEF_C,SPDEF_U), group = factor(rep(c("C","U"),each=5)))
leveneTest(c(CENPF_C,CENPF_U), group = factor(rep(c("C","U"),each=5)))
leveneTest(c(CDKN2A_C,CDKN2A_U), group = factor(rep(c("C","U"),each=5)))
leveneTest(c(E2F1_C,E2F1_U), group = factor(rep(c("C","U"),each=5)))
leveneTest(c(KLF2_C,KLF2_U), group = factor(rep(c("C","U"),each=5)))

rm(list = ls())
setwd("/data/nas1/lijia/58_KYGW-61101-6-NKT-KM112")

if (!dir.exists("06_正态性检验和方差齐性检验/")) {dir.create("06_正态性检验和方差齐性检验/")}
setwd("06_正态性检验和方差齐性检验/")


# ====================== 1. 录入数据 ======================
SPDEF_C <- c(0.942682752,0.954691221,1.132991853,1.050604248,0.919029926)
SPDEF_U <- c(1.182943271,1.227403386,1.292882879,1.07254862,1.335938059)

CENPF_C <- c(0.981486052,0.994218646,1.059037371,0.995606642,0.969651288)
CENPF_U <- c(1.000695423,1.057949833,1.30827587,1.107519324,1.124727116)

CDKN2A_C <- c(1.082723382,0.885031026,0.881823105,0.882340226,1.268082261)
CDKN2A_U <- c(1.177983854,1.003619284,0.936506515,1.039271396,1.363618962)

E2F1_C <- c(0.947360448,0.891899868,1.022410133,0.942612476,1.195717075)
E2F1_U <- c(1.070144674,1.260329761,1.083063529,1.211922882,1.297974981)

KLF2_C <- c(1.0830909,0.963966005,1.074347987,0.968165132,0.910429977)
KLF2_U <- c(0.692226314,0.607282024,0.770283685,0.68205064,0.779876053)

# 基因列表
genes <- c("SPDEF", "CENPF", "CDKN2A", "E2F1", "KLF2")

# 分组列表
groups <- list(
  list(C=SPDEF_C, U=SPDEF_U),
  list(C=CENPF_C, U=CENPF_U),
  list(C=CDKN2A_C, U=CDKN2A_U),
  list(C=E2F1_C, U=E2F1_U),
  list(C=KLF2_C, U=KLF2_U)
)

# ====================== 2. 批量做检验 ======================
library(car)

# 结果容器
res <- data.frame()

for (i in 1:5) {
  g <- genes[i]
  C <- groups[[i]]$C
  U <- groups[[i]]$U
  
  # 正态性检验
  sw_C <- shapiro.test(C)
  sw_U <- shapiro.test(U)
  
  # 方差齐性 Levene
  dt <- data.frame(value = c(C, U), group = rep(c("C","U"), each=5))
  lv <- leveneTest(value ~ group, data = dt)
  levene_p <- lv$`Pr(>F)`[1]
  
  # 结论
  normal_C <- ifelse(sw_C$p.value > 0.05, "正态", "非正态")
  normal_U <- ifelse(sw_U$p.value > 0.05, "正态", "非正态")
  levene_ok <- ifelse(levene_p > 0.05, "方差齐", "方差不齐")
  
  # 存入表格
  res <- rbind(res, data.frame(
    基因 = g,
    对照组_正态P = round(sw_C$p.value,4),
    对照组_分布 = normal_C,
    UCEC组_正态P = round(sw_U$p.value,4),
    UCEC组_分布 = normal_U,
    Levene方差齐P = round(levene_p,4),
    方差齐性 = levene_ok
  ))
}

# ====================== 3. 输出最终表格 ======================
cat("\n========== RT-qPCR 正态性与方差齐性检验结果表 ==========\n")
print(res, row.names = F)



# 保存模型比较结果到文本文件
sink("RESULTS.txt", append = FALSE)

cat("====================== RESULTS ======================\n\n")

cat(">init:\n")
print(res, row.names = F)
cat("\n")


cat("============================= END OF REPORT =============================\n")


# 结束保存
sink()