library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridBase)
library(RColorBrewer)
library(colorspace)

data_f <- snakemake@input[["cli"]]
cols_f <- snakemake@input[["colori"]]
cerchi <- snakemake@output[["circos"]]
leg <- snakemake@output[["legenda"]]

#dir <- '~/Dropbox/work/biobanca'
#egrassi@godot:/scratch/trcanmed/biobanca/local/share/data$ cat clinical_data_circos.tsv  | sed 's/TRUE/True/' | sed 's/FALSE/False/' > clinical_data_circos2.tsv

## for not removed cases use /scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/clinical_data_for_circos.tsv instead use
# /scratch/trcanmed/biobanca/local/share/data/complete_data_for_circos.tsv
#data <- read.table("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/clinical_data_for_circos.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
#data <- read.table("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/clinical_data_for_circos_revision.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
data <- read.table(data_f, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
for (i in rownames(data)) { 
  if (data[i, "buoni"] == "Validation failed") {
    data[i, "buoni"] <- "Failed"
  } else if (data[i, "buoni"] == "Validation not performed") {
    data[i, "buoni"] <- "BNot performed"
  } else if (data[i, "buoni"] == "Validation successful") {
    data[i, "buoni"] <- "Successful"
  } else {
    data[i, "buoni"] <- "BNot performed"
  }
}

names(data)[names(data)=="AGE.AT.COLLECTION..years."] <- "AGE.AT.COLLECTION"
names(data)[names(data)=="THERAPY.BEFORE.COLLECTION..Y.N."] <- "THERAPY.BEFORE..Y.N."

#data$Classification_N <- as.character(data$Classification_N)
#data$Classification_T <- as.character(data$Classification_T)
#cols <- read.table(file.path(dir, 'clinical_data_circos_cols_emendati_complete_difforder3_wderivation_revision.tsv'), sep="\t", header=TRUE, stringsAsFactors = FALSE, comment.char = "")
cols <- read.table(cols_f, sep="\t", header=TRUE, stringsAsFactors = FALSE, comment.char = "")
cols[is.na(cols)] <- "NA"


data[is.na(data$AGE.AT.COLLECTION),'AGE.AT.COLLECTION'] <- median(data$AGE.AT.COLLECTION, na.rm = TRUE)

outcircles <- c('circle.pdf')
outlegend <- c('circlelegend.pdf')

cols$col <- as.factor(cols$col)

create_named_vector <- function(level, data_cols) {
  my_col <- data_cols[data_cols$col == level,]
  if (unique(my_col$value)[1] == "numeric") { # we cannot have a value == numeric with this choice
    dir_palette <- strsplit(my_col$color, "_")
    res <- dir_palette # colorRamp2(breaks=c(30,50,60,70,80), colors=brewer.pal(n = 5, name = "RdBu"))
    names(res) <- 'numeric'
  } else if (unique(my_col$value)[1] == "automatic") { # neither automatic
    res <- my_col$color
    names(res) <- 'automatic'
  } else {
    res <- my_col$color
    names(res) <- my_col$value
  }
  return(res)
}

get_order <- function(level, data_cols) {
  my_col <- data_cols[data_cols$col == level,]
  res <-  unique(my_col$order)
  if (length(res) != 1) {
    stop('Each col in cols must be associated with a single order!')
  }
  return(res)
}


list_cols <- lapply(levels(cols$col), create_named_vector, cols)
names(list_cols) <- levels(cols$col)

order_cols <- sapply(levels(cols$col), get_order, cols)

########## Code to have a reasonable n. of samples
#tot <- data
#j <- 0
#for (i in seq(1, 25)) {
#  rownames(data) <- paste0(rownames(data), rep(j, nrow(data)))
#  j <- j+1
#  tot <- rbind(tot, data) 
#}
#data <- tot
#rownames(data) <- make.unique(replicate(nrow(data), paste(sample(LETTERS, 7), collapse="")))
#########

## error if we have annotation for a column not in data
stopifnot(length(levels(cols$col)) == length(intersect(colnames(data), levels(cols$col))))

circos.clear()

# We want to proceed in the user defined order to sort our complete data
order_cols <- order_cols[order(order_cols)]
#data <- data[order(data$sex, rownames(data)),]
order_columns <- paste(sapply(names(order_cols), function(x){paste0('data$', x)}), collapse=", ")
order_command <- paste0('data[order(', order_columns, ', rownames(data)),]')
data <- eval(parse(text=order_command))

# I am not sure that circlize is happy to work with *apply so we use a for (did not check though).
th <- (1/length(order_cols)) / 1.5
legends <- list()

pdf(cerchi)
for (i in seq(1, length(order_cols))) {
  name <- names(order_cols)[i]
  print(name)
  col_i <- list_cols[[name]]
  data_i <- data[, name, drop=FALSE]
  if (unique(names(col_i))[1] == "numeric") {
    # We need to build something like this:
    # colorRamp2(breaks=c(30,50,60,70,80), colors=brewer.pal(n = 5, name = "RdBu"))
    minv <- min(data_i[,1])
    maxv <- max(data_i[,1])
    q1 <- quantile(data_i[,1], 0.25, na.rm=TRUE)
    q3 <- quantile(data_i[,1], 0.75, na.rm=TRUE)
    medianv <- median(data_i[,1], na.rm=TRUE)
    #brewer <- brewer.pal(n=5, name=col_i[[1]][2])
    brewer <- sequential_hcl(5, col_i[[1]][2])
    if (col_i[[1]][1] == "minus") { # we may want to reverse the palette in some cases
      print(brewer)
      brewer <- rev(brewer)
      print(brewer)
    } 
    col_i <- colorRamp2(breaks=c(minv, q1, medianv, q3, maxv), colors=brewer)
    legends <- c(legends, Legend(title=name, col_fun=col_i, title_gp=gpar(fontsize=10, fontface="bold"), border = "black"))
  } else if (unique(names(col_i))[1] == "automatic") {
    my_factor <- as.factor(data_i[,1])
    my_levels <- levels(my_factor)
    palette <- col_i[[1]]
    col_i <- brewer.pal(n=length(my_levels), name=palette)
    names(col_i) <- my_levels
    legends <- c(legends, Legend(title=name, at=my_levels, legend_gp=gpar(fill=col_i), title_gp=gpar(fontsize=10, fontface="bold"), border = "black"))
  } else {
    my_factor <- as.factor(data_i[,1])
    my_levels_annot <- levels(my_factor)
    my_levels_conf <- levels(as.factor(names(list_cols[[name]])))
    delta <- setdiff(my_levels_annot, my_levels_conf)
    if (length(delta) != 0) {
      stop(paste0('I have values in column ', name, ' that have no annotation: ', delta, '\n'))
    }
    legends <- c(legends, Legend(title=name, at=names(list_cols[[name]]), legend_gp=gpar(fill=col_i), title_gp=gpar(fontsize=10, fontface="bold"), border = "black"))
  }
  #if (i == 1) {
  #  circos.heatmap(data_i, bg.border="black", bg.lwd=0.7, col=col_i, track.height=th, rownames.side="outside", rownames.cex=0.45)
  # } else{
  circos.heatmap(data_i, bg.border="black", bg.lwd=0.7, col=col_i, track.height=th)
  #}
}
graphics.off()

all_legends <- packLegend(list=legends, direction="horizontal", gap=unit(0.5, "cm"))
pdf(leg, width = 15, height = 10)
draw(all_legends)
graphics.off()
