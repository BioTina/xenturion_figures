#!/usr/bin/env Rscript
set.seed(42)

library(tidyverse)
library(networkD3)
library(htmlwidgets)


input <- snakemake@input[[1]]
#output_pdf <- snakemake@output[['pdf']]
output_html <- snakemake@output[['html']]

save.image("pippo.Rdata")

### the script takes one input: a df with ID-classification1-classification2
df <- read.table(input, sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
df[is.na(df)] <- "NC" 
df[,1] <- as.factor(df[,1])
df[,2] <- as.factor(df[,2])

links <-
    df %>% 
    mutate(row = row_number()) %>%  # add a row id
    gather('col', 'source', -row) %>%  # gather all columns
    mutate(col = match(col, names(df))) %>%  # convert col names to col nums
    mutate(source = paste0(source, '_', col)) %>%  # add col num to node names
    group_by(row) %>%
    arrange(col) %>%
    mutate(target = lead(source)) %>%  # get target from following node in row
    ungroup() %>% 
    filter(!is.na(target)) %>%  # remove links from last column in original data
    select(source, target) %>% 
    group_by(source, target) %>% 
    summarise(value = n())

# create nodes data frame from unque nodes found in links data frame
nodes <- data.frame(id = unique(c(links$source, links$target)),
                    stringsAsFactors = TRUE)

# remove column id from names
nodes$name <- as.factor(sub('_[0-9]*$', '', nodes$id))

# set links data to the 0-based index of the nodes in the nodes data frame
links$source <- match(links$source, nodes$id) - 1
links$target <- match(links$target, nodes$id) - 1

sn <- sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
              Target = 'target', Value = 'value', NodeID = 'name',
              fontSize = 12, nodeWidth = 30, iterations = 0)

saveWidget(sn, output_html)

library(webshot)
#webshot("sankey_model_cms.html", "sankey_model_cms.png")
#webshot(output_html, output_pdf)
