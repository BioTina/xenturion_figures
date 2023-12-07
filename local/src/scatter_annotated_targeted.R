#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
x <- args[1]
y <- args[2]
output <- args[3]
prefix <- args[4]
merge_col <- args[5]
yname <- args[6]
mut_f <- args[7]
cn_f <- args[8]
Rimage <- args[9]
output2 <- args[10]
annot_f <- args[11]
annot_w_f <- args[12]

load(Rimage)

data_x <- read.table(x, sep="\t", header=TRUE)
# y should be a two column file
data_y <- read.table(y, sep="\t", header=TRUE)
if (ncol(data_y) != 2 ) {
    stop(paste0(y, ' needs to be a two columns file!'))
}
colnames(data_y)[2] <- yname
data_x <- data_x[, c('case','CTG_5000')]
mdata <- merge(data_x, data_y, by=merge_col)
rownames(mdata) <- mdata[, merge_col]
mdata[,merge_col] <- NULL

cn <- read.table(cn_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
load(mut_f) # pdobing is a T/F matrix with genes on rows and models on cols TODO use only KNOWN MUTS! TODO
# here we have CRC0177 with the genotypes kept for molecular analyses (PDO is EGFR mutated, PDX is wt)
uniqued_mut <- t(pdobing)
uniqued_mut <- uniqued_mut[, colnames(uniqued_mut) %in% c('KRAS', 'BRAF', 'NRAS')]
uniqued_mut <- ifelse(uniqued_mut, 'MUT', 'WT')

#threewt <- as.data.frame(apply(uniqued_mut[,c('KRAS', 'BRAF', 'NRAS')], 1, function(x) {any(x!="wt")}))
#colnames(threewt) <- 'tf'
#threewt$smodel <- row.names(threewt)
#threewt$triple_wt <- ifelse(threewt$tf, 'MUT', 'WT')
#threewt$tf <- NULL
# merge_pdata2$triple_wt <- ifelse(infofour, 'triple MUT', 'WT')
#uniqued_mut2 <- as.data.frame(apply(uniqued_mut, 2, function(x) { ifelse(is.na(x),'wt', x) } ))
## there is a NA for PI3KCA but it has a KRAS mutation so that's fine for fourwt info

bin_cn <- as.data.frame(apply(cn, 2, function(x) {ifelse(x, 'AMPL', 'WT')}))
bin_cn <- bin_cn[, colnames(bin_cn) %in% c('ERBB2'), drop=FALSE]

save.image('pluto.Rdata')
uniqued_mut_w <- as.data.frame(uniqued_mut)
infofour <- apply(uniqued_mut_w[,c('KRAS', 'BRAF', 'NRAS')], 1, function(x) {any(x!="WT")})
infofour <- ifelse(is.na(infofour),  FALSE, infofour)
uniqued_mut_w$triple_wt <- ifelse(infofour, 'MUT', 'WT') 
dim(uniqued_mut_w)
uniqued_mut_w <- merge(uniqued_mut_w, bin_cn, by="row.names")
rownames(uniqued_mut_w) <- uniqued_mut_w$Row.names
uniqued_mut_w$Row.names <- NULL
dim(uniqued_mut_w)
uniqued_mut_w$alterations <- ifelse(uniqued_mut_w$ERBB2 != "WT", "ERBB2", uniqued_mut_w$triple_wt)
if (! all(uniqued_mut_w[uniqued_mut_w$ERBB2 != "WT", "triple_wt"] == "WT") ) {
  print('I cannot use a single color for ERBB2 and triple WT with whole!')
}
uniqued_mut_w[uniqued_mut_w$alterations == "MUT", 'alterations'] <- "KRAS/NRAS/BRAF"
uniqued_mut_w$smodel <- rownames(uniqued_mut_w)
write.table(uniqued_mut_w, file=annot_w_f, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

merge_pdata <- merge(mdata, bin_cn, by="row.names", all.x=TRUE)
rownames(merge_pdata) <- merge_pdata$Row.names
merge_pdata$Row.names <- NULL
merge_pdata2 <- merge(merge_pdata, uniqued_mut, by="row.names", all.x=TRUE)
# apply(merge_pdata2, 2, class) ?? character why
#infofour <- apply(merge_pdata2[,c('KRAS', 'NRAS', 'BRAF', 'PIK3CA')], 1, function(x) {any(x!="wt")})

# 3wt for Fra but not in the targeted, will in this way?
if (any(merge_pdata2$Row.names =="CRC1921")) { #manca nello shallow ma è 4wt per Fra
  merge_pdata2[merge_pdata2$Row.names=="CRC1921", 'NRAS'] <-  'WT'
  merge_pdata2[merge_pdata2$Row.names=="CRC1921", 'KRAS'] <-  'WT'
  merge_pdata2[merge_pdata2$Row.names=="CRC1921", 'BRAF'] <-  'WT'
}


infofour <- apply(merge_pdata2[,c('KRAS', 'BRAF', 'NRAS')], 1, function(x) {any(x!="WT")})
merge_pdata2$triple_wt <- ifelse(infofour, 'MUT', 'WT') 

if (! all(merge_pdata2[merge_pdata2$ERBB2 != "WT", "triple_wt"] == "WT") ) {
  stop('I cannot use a single color for ERBB2 and triple WT!')
}

merge_pdata2$alterations <- ifelse(merge_pdata2$ERBB2 != "WT", "ERBB2", merge_pdata2$triple_wt)
merge_pdata2[merge_pdata2$alterations == "MUT", 'alterations'] <- "KRAS/NRAS/BRAF"

write.table(merge_pdata2, file=annot_f, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
save.image('pippa.Rdata')
# there is a NA for PI3KCA but it has a KRAS mutation so that's fine
lmp <- function(model) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- model$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


ggplotRegression <- function(fit, data, name) {
  print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
          geom_point(data=data, aes(color=alterations), size=0.1) +
          stat_smooth(method = "lm", col = "darkgrey", size=0.2) +
          labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                             " P =",signif(summary(fit)$coef[2,4], 5)))+
          #theme_bw()+theme(text = element_text(size=15))+
          #scale_fill_manual(values=c('red', 'blue'))+
          #scale_colour_manual(values=c('red', 'black'))+
          scale_colour_manual(values=c('darkorange', 'darkred', 'grey40'))+
          #scale_shape_manual(values=c(15, 19))+
          #scale_size_manual(values=c(4, 2)) +
          #scale_x_continuous(breaks=c(0, 0.3, 0.6, 0.9, 1.2), limits=c(0, 1.2))+
          #scale_x_continuous(breaks=c(0.1, 0.3, 0.6, 0.9, 1.2))+
          labs(color="Relevant somatic alterations") + unmute_theme + # )
          theme(legend.position = "none"))
          #guides(fill = guide_legend(override.aes = list(shape = 19)),
          #       colour = guide_legend(override.aes = list(shape = 19))))
  #https://stackoverflow.com/questions/48043453/wrong-fill-values-in-a-ggplot2-legend
          #ggsave(paste0(name, ".svg"))
          ggsave(output, height=2.5, width=2.5, units='in')
}



compute_lm <- function(xcol, data, ycol) {
    model <- lm(data=data, formula=as.formula(paste0(ycol, "~", xcol)))
    ggplotRegression(model, data, xcol)
    cc <- cor.test(data[, xcol], data[,ycol])
    sm <- summary(model)
    r2 <- sm$r.squared
    ar2 <- sm$adj.r.squared
    pval <- lmp(sm) 
    ll <- logLik(model)
    return(c(pval, r2, ar2,ll, cc$estimate))
}

#compute_lm('CTG_5000', merge_pdata2, 'Cetuximab_dVw3')
res <- compute_lm('CTG_5000', merge_pdata2, yname)
#res <- t(sapply(xcolumns, compute_lm, mdata, yname))
res <- t(as.data.frame(res))
colnames(res) <- c('pvalue','R2','adjR2','logLik', 'pearson')
write.table(data.frame('exp'=rownames(res), res), file=output2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
