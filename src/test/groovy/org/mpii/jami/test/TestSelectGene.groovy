package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Shared
import spock.lang.Specification

/**
 * Created by mlist on 10/26/17.
 */
class TestSelectGene extends Specification {

    @Shared
    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")

    @Shared
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")

    @Shared
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")

    @Shared
    File testDir

    def setupSpec(){
        testDir = new File("out/test/")
        if (!testDir.exists()) {
            try {
                testDir.mkdirs()
            } catch (SecurityException se) {
                System.err.println("Could not create test directory.")
                System.err.println(se.getMessage())
                testDir = File.createTempDir()
            }
        } else if (!testDir.isDirectory()) {
            System.err.println("Could not create test directory.")
            testDir = File.createTempDir()
        }
    }

    def "test select gene"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_select_gene_filter.txt")

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
        def outputFileName = new File(testDir.absolutePath + "/test_select_gene_triplet.txt")

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
        def outputFileName = new File(testDir.absolutePath + "/test_select_gene_triplet.txt")

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
        def outputFileName = new File(testDir.absolutePath + "/test_select_gene_set.csv")
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