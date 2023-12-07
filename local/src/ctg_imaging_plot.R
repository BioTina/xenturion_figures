## scatter plot ctg imaging 
library(tidyverse)

dati_f <- snakemake@input[["pdo_cetuxi"]]
Rimage_f <- snakemake@input[["Rimage"]]
res_legend <- snakemake@output[["pic_legend"]]
res <- snakemake@output[["pic_square"]]
cors_f <- snakemake@output[["cors"]]

load(Rimage_f)
#dati_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/pdo_cetuxi_buoni.tsv"
dati <- read.table(dati_f, quote = "", sep = "\t", header = T, stringsAsFactors = F)
cor_5000 <- cor.test(dati$CTG_5000, dati$imaging_5000)
cor_1250 <- cor.test(dati$CTG_1250, dati$imaging_1250)
cor_20000 <- cor.test(dati$CTG_20000, dati$imaging_20000)

results <- as.data.frame(matrix(ncol = 2, nrow = 3))
rownames(results) <- c("cor_1250", "cor_5000", "cor_20000")
colnames(results) <- c("estimate_cor", "pvalue")
results[1,1] <- cor_1250$estimate
results[1,2] <- cor_1250$p.value
results[2,1] <- cor_5000$estimate
results[2,2] <- cor_5000$p.value
results[3,1] <- cor_20000$estimate
results[3,2] <- cor_20000$p.value

write.table(results, file=cors_f, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

gplot <- ggplot()+
  geom_point(data=dati, aes(x=CTG_1250, y=imaging_1250, colour = "1250"), size=0.1) + geom_smooth(data=dati, mapping = aes(x=CTG_1250, y=imaging_1250, colour = "1250", fill = "1250"), size=0.2, method = "lm", alpha = 0.1)+
  geom_point(data=dati, aes(x=CTG_5000, y=imaging_5000, color = "5000"), size=0.1) + geom_smooth(data=dati, mapping = aes(x=CTG_5000, y=imaging_5000, color = "5000", fill = "5000"), size=0.2, method = "lm", alpha = 0.1)+
  geom_point(data=dati, aes(x=CTG_20000, y=imaging_20000, color = "20000"), size=0.1) + geom_smooth(data=dati, mapping = aes(x=CTG_20000, y=imaging_20000, color = "20000", fill = "20000"), size=0.2, method = "lm", alpha = 0.1)+
  scale_color_manual(name = "Legend",
                     breaks = c("1250", "5000", "20000"),
                     values = c("1250" = "blue", "5000" = "orange", "20000" = "green"))+
  scale_fill_manual(name = "Legend",
                     breaks = c("1250", "5000", "20000"),
                     values = c("1250" = "blue", "5000" = "orange", "20000" = "green"))+
  xlab("CTG") + ylab("Imaging") + unmute_theme

ggsave(res_legend, plot=gplot, height=2.5, width=2.5, units='in')

gplot_nol <- gplot + theme(legend.position = "none")

ggsave(res, plot=gplot_nol, height=2.5, width=2.5, units='in')