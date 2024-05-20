library(ggplot2)
library(proxy)
library(pheatmap)

pheat_f <- snakemake@output[['pheat']]
#density_f <- snakemake@output[['density']]
#violin_f <- snakemake@output[['violin']]
#violin2_f <- snakemake@output[['violin2']]
jac_f <- snakemake@output[['jac_f']]
mw_f <- snakemake@output[['mw']]


load(snakemake@input[['Rimage']])

# we need to fill with 0 missing muts
thr <- 0.05
m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0
m1[m1 < thr] <- 0

mx <- m1[,grepl('.x', colnames(m1))]
mo <- m1[,grepl('.y', colnames(m1))]
colnames(mx) <- substr(colnames(mx), 0, 7)
colnames(mo) <- substr(colnames(mo), 0, 7)

# this removes CRC1367:
xx <- apply(mx, 2, sum)
oo <- apply(mo, 2, sum)
remove <- intersect(names(xx[xx==0]), names(oo[oo==0]))
print('Here we should have only CRC1367 remove for no muts, otherwise SOMETHING BAD HAS HAPPENED')
print(remove)
mx <- mx[, !colnames(mx) %in% remove]
mo <- mo[, !colnames(mo) %in% remove]


jac <- proxy::simil(mx, mo, by_rows = FALSE, method = "Jaccard")

pdf(pheat_f, family="sans")#, height=1.9, width=1.9)
ph <- pheatmap(jac, cluster_rows=F , cluster_cols=F, labels_col="PDXTs", labels_row="PDXs", fontsize.number=1.5, color=colorRampPalette(c("white", "red"))(50), legend=FALSE, border_color = "gray90" )
ph
graphics.off()
#pdf('urffa.pdf')
# found were to tweak it setting the color as darkgoldenrod and
# looking inside the objects
#ph$gtable[[1]][[1]]$children[[1]]$gp$lwd <- 0.001
#ph
#graphics.off()

write.table(jac, jac_f, sep="\t", quote=FALSE)
pearson <- jac
diag <- diag(pearson)
pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
sink(mw_f)
print('mean:')
mean(diag)
mean(all)
print('sd:')
sd(diag)
sd(all)
print('min:')
min(diag)
min(all)
print('max:')
max(diag)
max(all)
print('median:')
median(diag)
median(all)
print('summary:')
summary(diag)
summary(all)
sink()

pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))

#ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+geom_density(size=1.5)+scale_color_manual(values=c('#004D40','#FFC107'))+unmute_theme+xlab('jaccard')
#ggsave(density_f)


#ggplot(data=pdata[pdata$type=="matched",], aes(y=pearson,x=""))+geom_violin()+geom_jitter(height = 0, width = 0.1)+ylim(-0.001,1.0001)+unmute_theme+ylab('Mutational Similarity?')
#ggsave(violin_f)

#ggplot(data=pdata, aes(y=pearson,x=type, fill=type))+geom_violin()+geom_jitter(height = 0, width = 0.1)+ylim(-0.001,1.0001)+scale_fill_manual(values=c('#004D40','#FFC107'))+unmute_theme+xlab('jaccard')
#ggsave(violin2_f)

sink(mw_f, append=TRUE)
wilcox.test(formula=as.formula("pearson~type"), data=pdata)$p.value
sink()