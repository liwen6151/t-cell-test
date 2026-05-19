rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if (!dir.exists("07_Clinical")) {dir.create("07_Clinical")}
setwd("07_Clinical")

risk <- read.csv('../05_Prognostic_model/train/risk.csv', row.names = 1, check.names = F)
risk$id <- rownames(risk)
risk$id <- gsub('.', '-', risk$id , fixed = T)

train_dat <- read.csv('../00_rawdata/TCGA/01.fpkmlog2_TCGA.UCEC_mRNA.csv', row.names = 1, check.names = F)
max(train_dat)
train_dat <- train_dat[, colnames(train_dat) %in%risk$id]
colnames(train_dat) <- gsub('.', '-', colnames(train_dat), fixed = T)

survival <- read.csv("../00_rawdata/TCGA/01.survival_TCGA.UCEC.csv", row.names = 1, check.names = F)
survival$sample <- rownames(survival)
survival <- survival[survival$sample %in% colnames(train_dat), ]
survival$Sample <- substr(rownames(survival), 1, 12)

#######整理特征文件---------------------------------------
phenotype <- read.csv(file = "../00_rawdata/TCGA/data_clinical.csv")
colnames(phenotype)

train_phenotype <- data.frame(Sample = phenotype$bcr_patient_barcode,
                              age = phenotype$age_at_initial_pathologic_diagnosis,
                              grade = phenotype$neoplasm_histologic_grade,
                              #tumor_size = phenotype$intermediate_dimension.samples,
                              #TNM.stage=phenotype$tumor_stage.diagnoses
                              stage = phenotype$stage_event_clinical_stage,
                              # histologic = phenotype$histological_type,
                              tumor_invasion = phenotype$pct_tumor_invasion
)
train_phenotype <- merge(survival, train_phenotype, by = 'Sample')
train_phenotype <- subset(train_phenotype, select = -Sample)     
train_phenotype <- train_phenotype[, c("sample", setdiff(names(train_phenotype), "sample"))]
write.csv(train_phenotype, file = 'phenotype.csv')

train_phenotype$OS <- as.numeric(train_phenotype$OS)
train_phenotype$OS.time <- as.numeric(train_phenotype$OS.time)
train_phenotype2 <- train_phenotype


table(train_phenotype2$stage)
#train_phenotype2$stage <- gsub('stageajcc: ',' ',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IVB','IV',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IVA','IV',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IV','IV',train_phenotype2$stage, fixed = T)
#train_phenotype2$stage <- gsub('4','3',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIIB','III',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIIA','III',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIIC1','III',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIIC2','III',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIIC','III',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage III','III',train_phenotype2$stage, fixed = T)

#train_phenotype2$stage <- gsub('3','3_4',train_phenotype2$stage, fixed = T)

train_phenotype2$stage <- gsub('Stage IIB','II',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IIA','II',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage II','II',train_phenotype2$stage, fixed = T)
# train_phenotype2$stage <- gsub('Stage IB','1',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IA','I',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IB','I',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage IC','I',train_phenotype2$stage, fixed = T)
train_phenotype2$stage <- gsub('Stage I','I',train_phenotype2$stage, fixed = T)

#train_phenotype2$stage <- gsub('1','2',train_phenotype2$stage, fixed = T)
#train_phenotype2$stage <- gsub('2','1_2',train_phenotype2$stage, fixed = T)
train_phenotype2 <- subset(train_phenotype2, stage != "not reported")
table(train_phenotype2$stage)


table(train_phenotype2$grade)


table(train_phenotype2$tumor_invasion)
##空变成NA
# train_phenotype2$tumor_invasion <- ifelse(train_phenotype2$tumor_invasion == "", NA, train_phenotype2$tumor_invasion)
train_phenotype2$tumor_invasion <- gsub('tumor_invasion: ','',train_phenotype2$tumor_invasion, fixed = T)
train_phenotype2$tumor_invasion <- as.numeric(train_phenotype2$tumor_invasion)
train_phenotype2$tumor_invasion <- ifelse(train_phenotype2$tumor_invasion>=50,">=50","<50")
table(train_phenotype2$tumor_invasion)


table(train_phenotype2$age)
train_phenotype2$age <- gsub('age: ','',train_phenotype2$age, fixed = T)
train_phenotype2$age <- as.numeric(train_phenotype2$age)
train_phenotype2$age <- ifelse(train_phenotype2$age>=60,">=60","<60")
table(train_phenotype2$age)


train_phenotype3 <- train_phenotype2
train_phenotype2 <- train_phenotype3
colnames(train_phenotype2)



colnames(train_phenotype2) <- c('id', 'OS','OS.time',  'Age', 'Grade', 'Stage', "Tumor_invasion")
risk <- read.csv('../05_Prognostic_model/train/risk.csv', row.names = 1) %>% lc.tableToNum()
risk$id <- rownames(risk)
sub_risk <- subset(risk, select = c(id, riskScore,risk))
sub_risk$id <- gsub('.', '-', sub_risk$id , fixed = T)
train_risk_clinical <- merge(train_phenotype2, sub_risk, by = "id")

risk_clin <-train_risk_clinical
write.table(risk_clin,file="riskScore_Clincal.csv",sep=",",quote=F,row.names=F)





#画图--------------------
# risk_clin = train_risk_clinical
risk_clin <- na.omit(risk_clin)

risk_clin$riskScore <- as.numeric(risk_clin$riskScore)
#library(cowplot)
library(ggpubr)
sz=15
####riskScore####
########AGE
risk_clin$Age<-factor(risk_clin$Age,levels=c(">=60","<60"))
my_comparisons<-list(c(">=60","<60"))
g1<-ggboxplot(risk_clin, x="Age", y="riskScore", fill = "Age",
              ylab="riskScore",  xlab="Age",size = 1,axis.line =2)+
  theme(axis.text.y=element_text(size=sz),
        axis.title=element_text(size=sz),
        axis.text.x=element_text(size=sz)) +
  theme(plot.margin=unit(rep(3,4),'lines'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  stat_compare_means(comparisons=my_comparisons,label="p.signif",
                     method = "wilcox.test")
g1
pdf('riskScore_Clincal_Age.pdf')
print(g1)
dev.off()
png('riskScore_Clincal_Age.png')
print(g1)
dev.off()

#########PATHOLOGIC_Grade
table(risk_clin$Grade)

risk_clin_T<- risk_clin[risk_clin$Grade%in%c("G1","G2","G3","High Grade"),]
my_comparisons<-list(c("G1","G2"),c("G1","G3"),c("G1","High Grade"),c("G2","G3"),c("G2","High Grade"),c("G3","High Grade"))
risk_clin_T$Grade <- factor(risk_clin_T$Grade,levels =c("G1","G2","G3","High Grade"))
g3<-ggboxplot(risk_clin_T, x="Grade", y="riskScore", fill = "Grade",
              ylab="riskScore",  xlab="Grade",size = 1,axis.line =2)+
  theme(axis.text.y=element_text(size=sz),
        axis.title=element_text(size=sz),
        axis.text.x=element_text(size=sz)) +
  theme(plot.margin=unit(rep(3,4),'lines'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  stat_compare_means(comparisons=my_comparisons,label="p.signif")
# rotate_x_text(60)
g3
pdf('riskScore_Clincal_Grade.pdf')
print(g3)
dev.off()
png('riskScore_Clincal_Grade.png')
print(g3)
dev.off()

#########PATHOLOGIC_Stage
table(risk_clin$Stage)
risk_clin_N<- risk_clin[risk_clin$Stage%in%c("I","II","III","IV"),]
my_comparisons<-list(c("I","II"),c("I","III"),c("I","IV"),c("II","III"),c("II","IV"),c("III","IV"))
risk_clin_N$Stage <- factor(risk_clin_N$Stage,levels =c("I","II","III","IV"))

g4<-ggboxplot(risk_clin_N, x="Stage", y="riskScore", fill = "Stage",
              ylab="riskScore",  xlab="Stage",size = 1,axis.line =2)+
  theme(axis.text.y=element_text(size=sz),
        axis.title=element_text(size=sz),
        axis.text.x=element_text(size=sz)) +
  theme(plot.margin=unit(rep(3,4),'lines'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  stat_compare_means(comparisons=my_comparisons,label="p.signif")
g4
pdf('riskScore_Clincal_Stage.pdf')
print(g4)
dev.off()
png('riskScore_Clincal_Stage.png')
print(g4)
dev.off()

#########PATHOLOGIC_tumor_invasion
table(risk_clin$Tumor_invasion)
risk_clin_M<- risk_clin[risk_clin$Tumor_invasion%in%c("<50",">=50"),]
my_comparisons<-list(c("<50",">=50"))

g5<-ggboxplot(risk_clin_M, x="Tumor_invasion", y="riskScore", fill = "Tumor_invasion",
              ylab="riskScore",  xlab="Tumor_invasion",size = 1,axis.line =2)+
  theme(axis.text.y=element_text(size=sz),
        axis.title=element_text(size=sz),
        axis.text.x=element_text(size=sz)) +
  theme(plot.margin=unit(rep(3,4),'lines'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  stat_compare_means(comparisons=my_comparisons,label="p.signif")
g5
pdf('riskScore_Clincal_Tumor_invasion.pdf')
print(g5)
dev.off()
png('riskScore_Clincal_Tumor_invasion.png')
print(g5)
dev.off()
