library(dplyr)
library(ggplot2)
library(infotheo)
library(foreach)
library(iterators)

`10_genes_gene_expr` <- read.delim("/local/home/mlist/Projects/JAMI/data/10_genes_gene_expr.txt")
`10_genes_mir_expr` <- read.delim("/local/home/mlist/Projects/JAMI/data/10_genes_mir_expr.txt")
`10_genes_mirna_interactions_triplet_format` <- read.delim("/local/home/mlist/Projects/JAMI/data/10_genes_mirna_interactions_triplet_format.txt")

test_triplets_max_cores <- read.delim("/local/home/mlist/Projects/JAMI/out/test/test_triplets_max_cores.csv")
View(test_triplets_max_cores)

load("/local/home/mlist/Projects/JAMI/data/10_genes_cupid_results.RData")

comparison_data <-
    dplyr::left_join(test_triplets_max_cores, cupid_results,
                 by = c("Source" = "Modulator", "Target" = "Target", "miRNA" = "miRNA"),
                 suffix = c(".JAMI", ".CUPID"))


ggplot(comparison_data, aes(x = CMI.JAMI, y = CMI.CUPID)) + geom_point() +
    theme_bw() +
    geom_smooth(method=lm, se=FALSE) +
    annotate(geom = "text", x = 0.02, y = 0.1,
             label = paste("Pearson correlation:",
                       round(with(comparison_data,
                          cor(CMI.JAMI, CMI.CUPID, use = "complete.obs")), 2)))

ggplot(comparison_data, aes(x = p.value.JAMI, y = p.value.CUPID)) + geom_point() +
    theme_bw() +
    geom_smooth(method=lm, se=FALSE) +
    annotate(geom = "text", x = 0.1, y = 0.75,
             label = paste("Pearson correlation:",
                           round(with(comparison_data,
                                      cor(CMI.JAMI, CMI.CUPID, use = "complete.obs")), 2)))
