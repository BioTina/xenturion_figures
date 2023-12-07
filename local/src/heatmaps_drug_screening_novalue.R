library(reshape)
library(pheatmap)
library(readxl)

max_inh_f <- snakemake@input[['max_inh']]
targets_drug_screening <- snakemake@input[['drug_targets']]
heatmap_f <- snakemake@output[['heatmap']]

#d <- read.table('/scratch/trcanmed/biobanca/dataset/V1/drug_screening/max_drug_value.tsv', sep="\t", header=TRUE)
d <- read.table(max_inh_f, sep="\t", header=TRUE)
dm <- cast(d, formula=DRUG ~ MODEL, value = "MAX_VALUE", add.missing=TRUE, fill=NA)
row.names(dm) <- dm$DRUG
target <- data.frame(dm$DRUG)
dm$DRUG <- NULL

dm$mean <- rowMeans(dm, na.rm = TRUE)
dm$drug <- rownames(dm)
dm$drug2 <- sub("_[^.]*$", "", dm$drug)


targets_drug_screening <- read_excel("/scratch/trcanmed/biobanca/local/share/data/targets_drug_screening.xlsx", 
                                     col_names = FALSE)
colnames(targets_drug_screening) <- c("drug", "targets")
dm <- merge(dm, targets_drug_screening, by = "drug")
dm$name <- paste0(dm$drug2, "\n", dm$targets)
rownames(dm) <- dm$name
dm <- dm[order(-dm$mean),]

pheatmap(dm[2:6], cluster_cols=F, cluster_rows=F, file = heatmap_f)
