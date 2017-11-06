package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Specification

/**
 * Created by mlist on 10/26/17.
 */
class TestSelectGene extends Specification {

    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")

    def test_dir = new File("out/test").mkdir()

    def "test select gene"()
    {
        given:
        def outputFileName = new File("out/test/test_select_gene_filter.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = 100
        completeRun.filterForGene("ENSG00000110427")
        completeRun.runComputation()

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 40
    }

    def "test select gene triplets"()
    {
        given:
        def outputFileName = new File("out/test/test_select_gene_triplet.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = 100
        completeRun.selectedGenes = ["ENSG00000110427"]
        completeRun.runComputation()

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 40
    }

    def "test select two gene triplets"()
    {
        given:
        def outputFileName = new File("out/test/test_select_gene_triplet.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = 100
        completeRun.selectedGenes = ["ENSG00000110427", "ENSG00000100767"]
        completeRun.runComputation()

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 44
    }

    def "test select gene during reading set file"()
    {
        given:
        def outputFileName = new File("out/test/test_select_gene_set.csv")
        def genesMiRNA = new File("data/10_genes_mirna_interactions_set_format.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = 100
        completeRun.selectedGenes = ["ENSG00000110427"]
        completeRun.tripleFormat = false
        completeRun.runComputation()

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 40
    }
}