package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import spock.lang.Specification


/**
 * Created by mlist on 10/23/17.
 */
class TestParallelizationVsNoParallelizationVsNaive extends Specification {

    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")

    def test_dir = new File("out/test").mkdir()

    def "iterative partitioning with 1 core"() {
        given:
        def outputFileName = new File("out/test/test_triplets_1_core.csv")

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
        def outputFileName = new File("out/test/test_triplets_2_cores.csv")

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
        def outputFileName = new File("out/test/test_triplets_max_cores.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

    def "1 core and CUPID implementation"() {
        given:
        def outputFileName = new File("out/test/test_triplets_cupid.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        completeRun.method = "cupid";
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
    }

}