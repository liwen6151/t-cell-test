rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if(!dir.exists('./03_GO_KEGG_ppi/')){
  dir.create('./03_GO_KEGG_ppi')
}
setwd('./03_GO_KEGG_ppi/')

library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
options(stringsAsFactors = F)
library(ggplot2)
library(magrittr)
library(biomaRt)

#charset
go.circ = 1
kegg.circ = 1
kegg.buble = 1
kegg.chord = 1

# GO -------------------------------------------------------------
enrich_ID <- read.csv("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/02_hub_gene/01.gene_DGEs_TTKs_venn.csv", header = T)
select_DERNA <- read.csv("/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/01_DEGs/TCGA-UCEC/02.DEG_sig.csv",header=T)
if(str_detect(colnames(select_DERNA)[3],'logFC')){
  colnames(select_DERNA)[3] <- 'log2FoldChange'
}

gene_symbol <- bitr(
  geneID = enrich_ID$symbol,  #需要修改
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

enrich_go <- enrichGO(gene = gene_symbol[,2],   # 提供基因 ID 列
                      OrgDb = "org.Hs.eg.db",   # 使用人类基因注释数据库
                      keyType = "ENTREZID",     # 基因 ID 类型为 ENTREZID
                      ont = "ALL",              # 分析所有三种本体（BP、CC、MF）
                      pAdjustMethod = "fdr",    # 使用 FDR 方法调整 p 值
                      pvalueCutoff = 1,         # p 值截止阈值为 1
                      qvalueCutoff = 1,         # q 值截止阈值为 1
                      readable = TRUE)          # 返回可读的结果

enrich_go <- data.frame(enrich_go)
enrich_go <- mutate(enrich_go, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
GO_ALL <- enrich_go
if(median(GO_ALL$richFactor) < 0.1){GO_ALL$richFactor<- GO_ALL$richFactor*10}

GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]
paste0("得到",dim(GO_ALL)[[1]],"个结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
#"得到310个结果，其中262个生物学过程，29个细胞组分，19个分子功能"

write.csv(GO_ALL,file = "00.go_ALL.csv")
write.csv(as.data.frame(BP), '00.GO_BP.csv')
write.csv(as.data.frame(CC), '00.GO_CC.csv')
write.csv(as.data.frame(MF), '00.GO_MF.csv')

display_number = c(5, 5, 5)
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]
paste0(ego_result_BP$Description,collapse = "（）；")
paste0(ego_result_CC$Description,collapse = "（）；")
paste0(ego_result_MF$Description,collapse = "（）；")
col1 <- c("#E64B35B2", "#DC0000B2",'#E63863','#B2182B', '#D6604D','#EB4B17')
col2 <- c("#4DBBD5B2", "#0073C2B2",'#008280B2','#4393C3', '#23452F','#68A180')
col3 <- c("#E4C755", "#9FA3A8",'#D6E7A3','#F4A371FF','#F9C187FF','#FACF94FF')
Col1 <- sample(col1,1)
Col2 <- sample(col2,1)
Col3 <- sample(col3,1)

#GO柱状图-合并
if(T){
  go_enrich_df <- data.frame(
    ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
    Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
    type=factor(c(rep("BP:biological process", display_number[1]), 
                  rep("CC:cellular component", display_number[2]),
                  rep("MF:molecular function", display_number[3])), 
                levels=c("BP:biological process", "CC:cellular component","MF:molecular function" )))
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
  COLS <- c(Col1, Col2,Col3)
  
  p <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
    geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = COLS) + ###颜色
    coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "GO Terms")+
    theme_bw()+theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=12),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "italic", colour = "black")
    )+scale_x_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "01.go_barplot.pdf", height = 8, width = 12, p)
  ggsave(filename = "01.go_barplot.png", height = 8, width = 12, p)
} 

#circ
if(go.circ == 1){
  go_tem <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
  go <-  data.frame(Category = go_tem[,'ONTOLOGY'],ID = go_tem[,'ID'],Term = go_tem[,'Description'], 
                    Genes = gsub("/", ", ", go_tem[,'geneID']), adj_pval = go_tem[,'pvalue'])
  genelist <- data.frame(ID = select_DERNA$symbol, logFC = select_DERNA$log2FoldChange)
  row.names(genelist)=genelist[,1]
  library(GOplot)
  circ <- circle_dat(go, genelist)
  go_tem$pvalue <- -log10(go_tem$pvalue)
  go_enrich_df <- data.frame(
    category=go_tem$ONTOLOGY, 
    gene_num.min = 0,
    gene_num.max = max(go_tem$Count),
    gene_num.rich=go_tem$Count,
    "-log10.p" = go_tem$pvalue,
    up.regulated = 0,
    down.regulated = 0,
    rich.factor = go_tem$richFactor,
    Description=go_tem$ID)
  circ$sig <- ifelse(circ$logFC > 0,"up","down")
  tem1 <- go_enrich_df$Description
  i =1
  for(i in 1:15){
    x = tem1[i]
    circ1 <- circ[grep(x,circ$ID),]
    num1 <- nrow(circ1)
    tem2 <- as.data.frame(table(circ1$sig))
    if(nrow(tem2) == 2){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
    }else if(tem2[1,1] == "up"){
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]
    }else{go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]}
  }
  
  go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
  rownames(go_enrich_df) <- go_enrich_df$Description
  pdf(file = paste0("01.GOcirclize_plot.pdf"),width = 12,height = 6,family='Times')
  a <- dev.cur()   
  png(file = paste0("01.GOcirclize_plot.png"),width= 12, height= 6, units="in", res=300,family='Times')
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1,family="Times")
  circle_size = unit(1, 'snpc')
  #加载 circlize 包
  library(circlize)
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep(Col1, 5), rep(Col2, 5), rep(Col3, 5))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.6, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]
  label_data <- go_enrich_df['gene_num.rich'] 
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1 
  colorsChoice <- colorRampPalette(c('white','#4393C3'))
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  
  
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]
  label_data <- go_enrich_df['category']
  color_assign <- c('BP' = Col1, 'CC' = Col2, 'MF' = Col3)
  circos.genomicTrack(
    plot_data, ylim = c(0, max(plot_data$rich.factor)), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3) 
    } )
  circos.clear()
  category_legend <- Legend(
    labels = c("BP:biological process", "CC:cellular component","MF:molecular function" ),
    type = 'points', pch = NA, background = c(Col1,Col2,Col3), 
    labels_gp = gpar(fontsize = 12), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c('#FDDBC7', '#9ACD32'), 
    labels_gp = gpar(fontsize = 12), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 12), title_position = 'topleft', title = '-Log10(pvalue)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}


# KEGG--------------------------
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
options(stringsAsFactors = F)

gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$symbol,  
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

KEGG <- enrichKEGG( gene = gene_symbol$ENTREZID,#基因列表 
                    organism = "hsa",  #物种
                    keyType = "kegg",  #
                    minGSSize = 1, 
                    maxGSSize = 500,
                    pvalueCutoff = 1,  
                    pAdjustMethod = "BH",
                    qvalueCutoff = 1
)

kk <- setReadable(KEGG, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
KEGG <- data.frame(kk)
KEGG <- mutate(KEGG, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

kegg_result <- KEGG
if(median(kegg_result$richFactor) < 0.1){kegg_result$richFactor<- kegg_result$richFactor*10}
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
dim(hh)
# 12 10
hh <- hh[order(hh$Count,decreasing = T),]
write.csv(hh, file = "00.KEGG.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh[c(1:10),] 
paste0(kk$Description,collapse = "（）；")
if(kegg.circ == 1){
  go <-  data.frame(Category = "KEGG",ID = kk[,'ID'],Term = kk[,'Description'], 
                    Genes = gsub("/", ", ", kk[,'geneID']), adj_pval = kk[,'pvalue'])
  genelist <- data.frame(ID = select_DERNA$symbol, logFC = select_DERNA$log2FoldChange)
  row.names(genelist)=genelist[,1]
  library(GOplot)
  circ <- circle_dat(go, genelist)
  kk$pvalue <- -log10(kk$pvalue)
  go_enrich_df <- data.frame(
    category="KEGG", 
    gene_num.min = 0,
    gene_num.max = max(kegg_result$Count),
    gene_num.rich=kk$Count,
    "-log10.p" = kk$pvalue,
    up.regulated = 0,
    down.regulated = 0,
    rich.factor = kk$richFactor,
    Description=kk$ID)
  circ$sig <- ifelse(circ$logFC > 0,"up","down")
  tem1 <- go_enrich_df$Description
  for(i in 1:10){
    x = tem1[i]
    circ1 <- circ[grep(x,circ$ID),]
    num1 <- nrow(circ1)
    tem2 <- as.data.frame(table(circ1$sig))
    if(nrow(tem2) == 2){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
    }else if(tem2[1,1] == "up"){
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]
    }else{
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]}
  }
  
  go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
  rownames(go_enrich_df) <- go_enrich_df$Description
  library(ggplot2)
  pdf(file = paste0("02.KEGG_cic.pdf"),width = 13,height = 6,family='Times')
  a <- dev.cur()   
  png(file = paste0("02.KEGG_cic.png"),width= 13, height= 6, units="in", res=300,family='Times')
  dev.control("enable")
  par(family="Times")
  circle_size = unit(1, 'snpc')
  library(circlize)
  kegg.color <- sample(col1,1)
  circos.par(gap.degree = 2, start.degree = 90)
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep(kegg.color, 10))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.6, major.tick.length = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('KEGG' = kegg.color)
  circos.genomicTrack(
    plot_data, ylim = c(0, max(plot_data$rich.factor)), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value,col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    })
  circos.clear()
  category_legend <- Legend(
    labels = c("KEGG Pathway"),
    type = 'points', pch = NA, background = c(kegg.color), 
    labels_gp = gpar(fontsize = 12), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background =c('#FDDBC7', '#9ACD32'), 
    labels_gp = gpar(fontsize = 12), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 12), title_position = 'topleft', title = '-Log10(pvalue)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}


save.image('enrich.Rdata')


#ppi分析-----------
# 清除环境变量，设置工作目录
rm(list = ls())
workpath = '/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/'
setwd(workpath)
if (!dir.exists("./03_GO_KEGG_ppi/")) {dir.create("./03_GO_KEGG_ppi/")}
setwd('./03_GO_KEGG_ppi/')

# 加载R包
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
library(export)
library(circlize)
library(grid)
library(graphics)
library(ComplexHeatmap)
select=dplyr::select

# 爬虫
get.string <- function(ids, score = 400, species = 9606){
  
  string_api_url <- "https://string-db.org/api"
  # Mapping identifiers------
  # 9606: 人类 (默认)
  # 10090: 小鼠
  # 10116: 大鼠
  output_format <- "tsv"
  method <- paste0("get_string_ids?species=", species, "&identifiers=")
  identifiers <- paste0(ids, collapse = '%0d')
  url <- paste0(c(string_api_url, output_format, method), collapse = '/')
  url <- paste0(url, identifiers)
  tmp <- tempfile()
  curl::curl_download(url, tmp)
  identifiers <- read.delim2(tmp) %>% pull(preferredName) %>% paste0(collapse = '%0d')
  unlink(tmp)
  
  
  # Getting the STRING network interactions-----
  # 150：最低置信度（low confidence）0.15
  # 400：中等置信度（medium confidence, 默认）0.4
  # 700：高置信度（high confidence）0.7
  # 900：最高置信度（highest confidence）0.9
  output_format <- "tsv"
  method <- paste0("network?species=", species, "&required_score=", score, "&identifiers=")
  url <- paste0(c(string_api_url, output_format, method),collapse = '/')
  url <- paste0(url, identifiers)
  curl::curl_download(url, 'string_interactions.tsv')
  
  
  # Getting STRING network image-------
  output_format <- "highres_image" # 可选 "image","svg"
  method <- paste0("network?species=", species, "&required_score=", score, "&identifiers=") 
  url <- paste0(c(string_api_url, output_format, method), collapse = '/')
  url <- paste0(url,identifiers)
  # curl::curl_download(url,'string_hires_image.png')
  
  print("finished")
  
}

gene <- read.csv("../02_hub_gene/01.gene_DGEs_TTKs_venn.csv")

# 下载ppi网络信息
options(timeout = 9999999)
get.string(gene$symbol, score = 150, species = 9606)


string_data <- read_tsv('string_interactions.tsv') %>% dplyr::select(c("preferredName_A", "preferredName_B", "score"))

# string_data <- read_tsv('string_interactions_short.tsv') %>%
#   dplyr::select(c("preferredName_A"="#node1", "preferredName_B" = "node2", "score" = "combined_score")) %>%
#   filter(score>=0.7)

p_fun <- function(string_data,select_dregee=1){
  
  # 步骤1：初始化过滤过程
  filtered_data <- string_data  # 先使用原始数据
  has_changes <- TRUE  # 设置初始条件以进入循环
  
  # 步骤2：迭代过滤，直到不再有变化
  while(has_changes) {
    # 获取当前数据中的所有节点
    all_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
    
    # 计算每个节点的度数（连接数）
    node_degrees <- as.data.frame(table(all_nodes))
    colnames(node_degrees) <- c("node", "degree")
    
    # 找出度数≥2的节点
    nodes_to_keep <- node_degrees$node[node_degrees$degree >= select_dregee]
    
    # 过滤数据，只保留度数≥2的节点之间的连接
    new_filtered_data <- filtered_data %>%
      filter(preferredName_A %in% nodes_to_keep & 
               preferredName_B %in% nodes_to_keep)
    
    # 检查是否有变化
    if(nrow(new_filtered_data) == nrow(filtered_data)) {
      has_changes <- FALSE  # 没有变化，退出循环
    } else {
      filtered_data <- new_filtered_data  # 更新数据，继续迭代
    }
  }
  
  # 步骤3：检查结果
  print(paste("原始数据行数:", nrow(string_data)))
  print(paste("过滤后数据行数:", nrow(filtered_data)))
  
  # 最终确认所有保留的节点都有度数≥2
  final_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
  final_degrees <- as.data.frame(table(final_nodes))
  colnames(final_degrees) <- c("node", "degree")
  
  # 检查是否有度数<2的节点
  low_degree_nodes <- final_degrees$node[final_degrees$degree < select_dregee]
  if(length(low_degree_nodes) > 0) {
    print(paste0("警告：仍有度数<",select_dregee,"的节点:"))
    print(low_degree_nodes)
  } else {
    print(paste0("成功：所有节点的度数都≥",select_dregee))
  }
  
  print(paste0('包含节点：', length(unique(c(filtered_data$preferredName_A,filtered_data$preferredName_B)))))
  print(paste0('包含边：', nrow(filtered_data)))
  
  library(ComplexHeatmap)
  library(circlize)
  # 准备数据
  chord_data <- data.frame(
    from = filtered_data$preferredName_A,
    to = filtered_data$preferredName_B,
    value = filtered_data$score
  )
  
  # 获取唯一的基因名称
  all_genes <- unique(c(chord_data$from, chord_data$to))
  
  # 计算每个基因的连接度（与其他节点的连接数量）
  gene_connectivity <- table(c(chord_data$from, chord_data$to))
  gene_connectivity <- gene_connectivity[all_genes]  # 确保顺序一致
  
  # 按照连接度从高到低排序基因
  sorted_genes <- names(sort(gene_connectivity, decreasing = TRUE))
  
  # 创建连接度的颜色渐变 - 使用ComplexHeatmap的colorRamp2函数
  connectivity_range <- range(gene_connectivity)
  connectivity_col_fun <- colorRamp2(
    breaks = c(connectivity_range[1], mean(connectivity_range), connectivity_range[2]),
    colors = c("blue", "green", "red")
  )
  sector_colors <- connectivity_col_fun(as.numeric(gene_connectivity))
  names(sector_colors) <- all_genes
  
  # 创建score的颜色渐变函数
  score_range <- range(chord_data$value)
  score_col_fun <- colorRamp2(
    breaks = c(score_range[1], mean(score_range), score_range[2]),
    colors = c("lightblue", "purple", "darkred")
  )
  
  par(mar = c(1, 1, 3, 1))  # 增加顶部边距
  
  # 清除任何现有的circos布局
  circos.clear()
  
  # 设置基于扇区数量的较小间隔
  n_sectors <- length(all_genes)
  gap_degree <- min(2, 360/(n_sectors * 3))  # 自适应间隔大小
  
  # 设置圆形布局参数 - 从12点钟方向开始
  circos.par(
    gap.degree = gap_degree,
    start.degree = 90,  # 90度是12点钟方向
    # 调整画布大小，留出更多空间给标签
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1.2, 1.2)  # 上方留出更多空间
  )
  
  # 创建弦图，使用order参数按连接度排序
  chordDiagram(
    chord_data,
    grid.col = sector_colors,  # 根据连接度的颜色
    col = score_col_fun,  # 连线根据score着色
    transparency = 0.5,
    link.lwd = 1.5,
    link.lty = 1,
    link.sort = TRUE,
    link.decreasing = TRUE,
    annotationTrack = "grid",  # 只显示网格线，不显示标签
    order = sorted_genes  # 按照连接度排序基因
  )
  
  # 添加更外围的标签
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(c(get.cell.meta.data("xlim"))),
        y = get.cell.meta.data("ylim")[1] + 1.2,  # 将标签放置得更靠外
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        # cex = 0.3
        cex = 0.4
      )
    },
    bg.border = NA
  )
  title("", line = -2)
  # 添加标题 - 放置得更高
  title(
    "Protein-Protein Interaction Network",
    line = 0.3,  # 增加此值会将标题放置得更高
    cex.main = 1.2,  # 稍微增大标题字体
    font.main = 2    # 粗体标题
  )
  
  # 使用ComplexHeatmap的Legend函数创建图例
  # 为连接度创建图例
  connectivity_legend <- Legend(
    title = "Degree",
    col_fun = connectivity_col_fun,
    at = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    labels = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  # 为score值创建图例
  score_legend <- Legend(
    title = "Interaction Score",
    col_fun = score_col_fun,
    at = round(seq(score_range[1], score_range[2], length = 5), 3),
    labels = round(seq(score_range[1], score_range[2], length = 5), 3),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  # 将两个图例打包在一起
  legend_list <- packLegend(connectivity_legend, score_legend, direction = "vertical")
  
  # 绘制图例
  pushViewport(viewport(x = 0.85, y = 0.15, width = 0.3, height = 0.3))
  draw(legend_list)
  popViewport()
  
  # 重置circlize参数
  circos.clear()
  
}

# [1] "原始数据行数: 493"
# [1] "过滤后数据行数: 493"
# [1] "成功：所有节点的度数都≥1"
# [1] "包含节点：83"
# [1] "包含边：493"

library(Cairo)
CairoPDF("03.ppi.pdf", width=9, height=8)
p_fun(string_data)
dev.off()

CairoPNG("03.ppi.png", width=9, height=6.5, units="in", res=300)
p_fun(string_data)
dev.off()

