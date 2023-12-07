library(ggplot2)
library(ggrastr)

pdodata <- snakemake@input[["pdo"]]
pdxdata <- snakemake@input[["pdx"]]
Rimage_f <- snakemake@input[["Rimage"]]
result <- snakemake@output[["pdf"]]
#odata <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv')
#pdata <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv')
odata <- read.table(pdodata)
pdata <- read.table(pdxdata)
m <- merge(odata, pdata, by="row.names")

nrow(m) # for the caption, need to put in a log file
ci <- cor.test(m$log2FoldChange.x, m$log2FoldChange.y) # idem
ci
ci$p.value
#dim(m)
#dim(odata)
#dim(pdata)
load(Rimage_f)

m$padj_x <- ifelse(m$padj.x < 10**(-10), 10**(-10), m$padj.x)
p <- ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
  rasterize(geom_point(size=0.1), dpi=300)+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme+theme(legend.position = "none")
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
ggsave(p, file=result, width=4, height=4, units="in")

#ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
#  geom_point(size=0.1)+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
#  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
#ggsave("~/biobanca_cor_ctx_legend.pdf", width=4, height=4, units="in")

# 
# m$padj_x <- ifelse(m$padj.x < 10**(-10), 10**(-10), m$padj.x)
# ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
#   geom_point()+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
#   scale_color_distiller(palette = "YlOrRd", direction=1)+theme_bw()