package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Shared
import spock.lang.Specification

/**
 * Created by mlist on 10/26/17.
 */
class TestPValueCutoff extends Specification {

    @Shared
    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")

    @Shared
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")

    @Shared
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")

    @Shared
    def numberOfPermutations = 100

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

    def "test p-value cutoff"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_p_value_cutoff.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.pValueCutoff = 0.05
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk < 100
    }
}