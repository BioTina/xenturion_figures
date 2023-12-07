###pheatmap for correlation
library(pheatmap)
c_s <- snakemake@input[["correlation_matrix"]]
pheatmap_correlation <- snakemake@output[["pheatmap"]]
correlation_simo <- read.table(c_s, quote = "", sep = "\t", header = TRUE, row.names = 1)

pdf(pheatmap_correlation)#, width=10, height=8, units="in", type="cairo", res=300)
pheatmap(correlation_simo, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

