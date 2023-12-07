library(precrec)
library(ggplot2)
input_auc_f <- snakemake@input[['input_auc']]
targeted_annot_f <- snakemake@input[['targeted_annot']]
auc_all_f <- snakemake@output[['auc_all']]
auc_wt_f <- snakemake@output[['auc_wt']]
sens_f <- snakemake@output[['sens_table']]
log_f <- snakemake@log[['log']]

data <- read.table(input_auc_f, sep="\t", header=TRUE)
#data <- read.table('/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/input_auc.tsv', sep="\t", header=TRUE)
data$label <- data$responder
precrec_obj <- evalmod(scores = -data$CTG_5000, labels = data$label)
pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))
graphics.off()

sink(log_f)
precrec_obj
sink()

get_sensitivity_specificity <- function(thr, df) {
  p <- nrow(df[df$labels == 1,])
  tp <- nrow(df[df$scores > thr & df$labels==1,])
  n <- nrow(df[df$labels == 0,])
  tn <- nrow(df[df$scores <= thr & df$labels==0,])
  fn <- nrow(df[df$scores <= thr & df$labels==1,])
  fp <- nrow(df[df$scores > thr & df$labels==0,])
  return(c(tp/p, tn/n, tp, tn, fp, fn, fp/(fp+tp)))
}

compute_thr <- function(scores, labels) {
  df <- data.frame(scores=scores, labels=labels)
  intervals <- cut(scores, breaks=10)
  max_i <- levels(intervals)
  max_i <- gsub('(','', max_i, fixed=TRUE)
  max_i <- gsub(']','', max_i, fixed=TRUE)
  upper_bs <- sapply(strsplit(max_i, ','), '[[', 2)
  upper_bs <- as.numeric(upper_bs)
  upper_bs <- upper_bs[order(-upper_bs)]
  res <- as.data.frame(t(sapply(upper_bs, get_sensitivity_specificity, df)))
  colnames(res) <- c('sensitivity', 'specificity', 'tp', 'tn', 'fp', 'fn', 'fdr')
  res$thr <- upper_bs
  return(res)
}

sens <- compute_thr(-data$CTG_5000, data$label)
write.table(sens, file=sens_f, sep="\t", quote=FALSE, row.names=FALSE)
#annot <- read.table('/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/targeted_annot.tsv', sep="\t", header=TRUE)
annot <- read.table(targeted_annot_f, sep="\t", header=TRUE)
wt <- annot[annot$alterations=='WT', 'Row.names']
data <- data[data$case %in% wt, ]
precrec_obj <- evalmod(scores = -data$CTG_5000, labels = data$label)
pdf(auc_wt_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))
graphics.off()
sink(log_f, append=TRUE)
precrec_obj
sink()
