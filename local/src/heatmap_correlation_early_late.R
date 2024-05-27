library(tidyverse)
library(pheatmap)

vsd_f <- snakemake@input[["vsd_file"]]
samples_f <- snakemake@input[["samples_file"]]
pdf <- snakemake@output[["res_pdf"]]
#res_tsv <- snakemake@output[["result"]]

#vsd_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/early_late/vsd.tsv.gz"
vsd <- read.table(vsd_f, sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)

#samples_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/early_late/samples_data"
samples <- read.table(samples_f, quote = "", sep = "\t",header = TRUE, stringsAsFactors = FALSE)

#vsd <- as.data.frame(t(vsd))
#vsd$id <- rownames(vsd)
#merged <- merge(vsd, samples, by="id")
#new <- merged
#new$id <- NULL

### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
means <- apply(vsd, 1, mean)
med <- median(means)

de <- vsd[means > med,]

sds <- apply(de, 1, sd)

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

desd1 <- desd
#names(desd1) <- substr(names(desd1), 1, 12)
desd1 <- as.data.frame(t(desd1))
desd1$id <- rownames(desd1)
merged <- merge(desd1, samples, by="id")
new <- merged
new$model_passage <- paste0(new$model, "_", new$passage)
new$id <- NULL
new$model <- NULL
new$passage <- NULL
rownames(new) <- new$model_passage
new$model_passage <- NULL
new <- as.data.frame(t(new))

early <- new[grepl('early', names(new), , fixed=TRUE)]
late <- new[grepl('late', names(new), , fixed=TRUE)]


colnames(early) <- substr(colnames(early),0,7)
colnames(late) <- substr(colnames(late),0,7)
cearly <- early[, names(early)]
clate <- late[, names(cearly)]
### check also for genes
all_genes <- intersect(rownames(early), rownames(late))
cearly <- cearly[all_genes,]
clate <- clate[all_genes,]


if (all(colnames(cearly)!=colnames(clate)) & all(rownames(cearly)!=rownames(clate))) {
  stop('Brutto llama!')
}

colnames(cearly) <- paste0(colnames(cearly), "_early")
colnames(clate) <- paste0(colnames(clate), "_late")

res <- cor(cearly, clate)

pdf(file=pdf)
pheatmap(res, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()
#write.table(res, file=res_tsv, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

diagonal_values <- diag(res)
#quantile(diagonal_values)

non_diagonal_values <- res[lower.tri(res) | upper.tri(res)]
#quantile(non_diagonal_values)

# print("diag")
# quantile(diagonal_values)
# print("non_diag")
# quantile(non_diagonal_values)