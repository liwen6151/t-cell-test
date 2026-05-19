#hub_gene_TTK-------------------------------------
rm(list = ls())
setwd('/data/nas1/zhangzhaolei/project/05.KYGW-61101-6-NKT/')
if (! dir.exists("./02_hub_gene/")){
  dir.create("./02_hub_gene")
}
setwd("./02_hub_gene")

gene_A <- read.csv('../01_DEGs/TCGA-UCEC/02.DEG_sig.csv')
gene_B <- read.csv('../00_rawdata/TTK_gene.csv')

inter_gene <- intersect(gene_A$symbol, gene_B$symbol) %>% as.data.frame()
colnames(inter_gene) <- 'symbol'
write.csv(inter_gene, '01.gene_DGEs_TTKs_venn.csv')


library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)
library(ggtext)

venn_list <- list(DEGs = gene_A$symbol,
                  TTK = gene_B$symbol)

txt <- data.frame(
  x = c(0),
  y = c(1.8),
  label =paste0('Complex Venn Diagram')) 
p <- ggvenn(venn_list, 
            c('DEGs', "TTK"),
            fill_color = c("#FFD306","#871F78"),
            show_percentage = T, # 显示交集数量的百分比
            stroke_alpha = 0.5,
            stroke_size = 0.5, # 交集处白边的大小
            stroke_color="white",
            stroke_linetype="solid",
            text_size = 5,
            set_name_color=c("#FFD306","#871F78"),
            text_color = 'black')
p
ggsave('02.hub_gene.pdf', p, w=8, h=8)
ggsave('02.hub_gene.png', p, w=8, h=8)

