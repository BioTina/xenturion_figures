### Draw densit plot
library(ggplot2)
library(pheatmap)
res_f<- snakemake@input[["res"]]
density_f <- snakemake@output[["density_plot"]]
wilcox_f <- snakemake@output[["wilcox_result"]]

#theme ggplot 

textSize <- 5
largerSize <- textSize + 2

unmute_theme <- theme_bw() +
theme(
	text = element_text(size = textSize, family='sans'),
	axis.title = element_text(size = largerSize),
	axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
	axis.text.y = element_text(size = textSize, color="black"),
	plot.title = element_text(size = largerSize, hjust = 0.5),
	legend.title = element_text(size=largerSize),
    legend.text = element_text(size=textSize)
)


res_df <- read.table(gzfile(res_f), quote = "", sep = "\t", header = TRUE)

### prendo la diagonale, quindi i campioni matched con la funzione diag mentre gli unmatched li prendo con upper.tri per la
### parte sopra la diagonale e lower.tri per la parte inferiore alla diagonale

matched <- diag(as.matrix((res_df))) ### diag vuole la matrice perchè non converte da solo il data frame
unmatched <- c(res_df[upper.tri(res_df)], res_df[lower.tri(res_df)])

res_data <- data.frame(res_df=c(unmatched, matched), type=c(rep('unmatched', length(unmatched)), rep('matched', length(matched))))
ggplot(data=res_data, aes(x=res_df, color=type))+geom_density()+xlim(c(0,1))+scale_color_manual(values=c("darkgreen", "darkgoldenrod1"))+unmute_theme+xlab("pearson")+theme(legend.position="none")
ggsave(density_f)

#write.table(as.data.frame(pearson), file=pearson_f, quote=FALSE, sep="\t")
#save.image(rdata_f)

### test wilkcoxson
x <- matched
y <- unmatched
w <- wilcox.test(x, y)

wilcox <- data.frame(row.names = "wilcox.test", pvalue=w$p.value)
write.table(wilcox, file=wilcox_f, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

sink(wilcox_f, append=TRUE)
length(x)
length(y)
summary(matched)
summary(unmatched)
sink()

