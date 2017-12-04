package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import org.mpii.jami.model.Triplet
import spock.lang.Shared
import spock.lang.Specification

import static spock.util.matcher.HamcrestMatchers.closeTo


/**
 * Created by mlist on 10/24/17.
 */
class TestVariousCmiImplementations extends Specification {

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

    @Shared
    Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")


    def "test uniform grid"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_uniform.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.numberOfBins = 3
        completeRun.method = "uniform"
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.038599, 0.001)
    }

    def "test pseudouniform grid"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_pseudouniform.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.numberOfBins = 3
        completeRun.method = "pseudouniform"
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.0809767, 0.001)
    }

    def "test cupid"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_cupid.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.method = "cupid"
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.08391, 0.001)
    }

    def "test normal"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_iterative_partitioning.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.14611, 0.001)
    }

    def "test normal with set format"()
    {
        given:
        genesMiRNA = new File("data/10_genes_mirna_interactions_set_format.txt")
        def outputFileName = new File(testDir.absolutePath + "/test_iterative_partitioning_set.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName)
        completeRun.numberOfPermutations = numberOfPermutations
        completeRun.tripleFormat = false
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 342
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.14611, 0.001)
    }

}