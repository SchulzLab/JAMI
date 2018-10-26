package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Shared
import spock.lang.Specification


/**
 * Created by mlist on 10/23/17.
 */
class TestRandomSeeding extends Specification {

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

    def "repeated iterative partitioning with 2 cores and no seed"() {
        given:
        def outputFileNameA = new File(testDir.absolutePath + "/test_triplets_no_seed_a.txt")
        def outputFileNameB = new File(testDir.absolutePath + "/test_triplets_no_seed_a.txt")

        when:
        CompleteRun completeRunA = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileNameA)
        CompleteRun completeRunB = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileNameB)

        completeRunA.numberOfThreads = 2
        completeRunB.numberOfThreads = 2

        completeRunA.runComputation()
        completeRunB.runComputation()

        then:
        completeRunA.completed == true
        completeRunB.completed == true

        completeRunA.triplets.collect{it.pValue} != completeRunB.triplets.collect{it.pValue}
    }

    def "repeated iterative partitioning with 2 cores and seed"() {
        given:
        def outputFileNameA = new File(testDir.absolutePath + "/test_triplets_seed_a.txt")
        def outputFileNameB = new File(testDir.absolutePath + "/test_triplets_seed_a.txt")

        when:
        CompleteRun completeRunA = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileNameA)
        CompleteRun completeRunB = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileNameB)

        completeRunA.numberOfThreads = 2
        completeRunB.numberOfThreads = 2

        completeRunA.setSeed(12345)
        completeRunB.setSeed(12345)

        completeRunA.runComputation()
        completeRunB.runComputation()

        then:
        completeRunA.completed == true
        completeRunB.completed == true

        completeRunA.triplets.collect{it.pValue} == completeRunB.triplets.collect{it.pValue}
    }
}