package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Shared
import spock.lang.Specification


/**
 * Created by mlist on 10/23/17.
 */
class TestParallelizationVsNoParallelizationVsNaive extends Specification {

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


    def "iterative partitioning with 1 triplet and 1 core"() {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_triplets_1_core_1_interaction.txt")
        def genesMiRNA_one = new File("data/10_genes_mirna_interactions_triplet_format_1_triplet.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA_one,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfThreads = 1;
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 4
    }

    def "iterative partitioning with 1 core"() {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_triplets_1_core.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfThreads = 1;
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

    def "iterative partitioning with 2 cores"() {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_triplets_2_cores.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        completeRun.numberOfThreads = 2;
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

    def "iterative Partitioning -1 cores"() {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_triplets_max_cores.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

    def "1 core and CUPID implementation"() {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_triplets_cupid.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        completeRun.method = "cupid";
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

}