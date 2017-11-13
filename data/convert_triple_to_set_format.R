library(foreach)
library(dplyr)

`125_genes_ceRNA_interactions_to_test_in_cupid` <- read.delim("/local/home/mlist/Projects/JAMI/data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
View(`125_genes_ceRNA_interactions_to_test_in_cupid`)

set_format <- foreach(gene = unique(c(
    as.character(`125_genes_ceRNA_interactions_to_test_in_cupid`$geneA),
    as.character(`125_genes_ceRNA_interactions_to_test_in_cupid`$geneB))),
    .combine = rbind) %do%
    {
        data.frame(gene = gene, miRNA = paste(
            unique(as.chCaracter(dplyr::filter(`125_genes_ceRNA_interactions_to_test_in_cupid`,
                      geneA == gene | geneB == gene)$miRNA)), collapse = ","))
    }

write.table(set_format, row.names = FALSE, quote = FALSE, sep = "\t",
            file = "~/Projects/JAMI/data/125_genes_ceRNA_interactions_set_format.txt")
