library(foreach)

`125_genes_ceRNA_interactions_to_test_in_cupid` <- read.delim("/local/home/mlist/Projects/JAMI/data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
View(`125_genes_ceRNA_interactions_to_test_in_cupid`)

foreach(gene = unique(c(`125_genes_ceRNA_interactions_to_test_in_cupid`$geneA,
                 `125_genes_ceRNA_interactions_to_test_in_cupid`$geneB)), .combine = rbind) %do%
    {
        data.frame(gene = gene, miRNA = paste(unique(dplyr::filter(`125_genes_ceRNA_interactions_to_test_in_cupid`,
                      geneA == gene | geneB == gene)$miRNA), collapse = ","))
    }
