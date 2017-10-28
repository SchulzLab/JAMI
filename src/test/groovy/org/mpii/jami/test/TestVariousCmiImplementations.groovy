package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import org.mpii.jami.model.Triplet
import spock.lang.Specification

import static spock.util.matcher.HamcrestMatchers.closeTo


/**
 * Created by mlist on 10/24/17.
 */
class TestVariousCmiImplementations extends Specification {

    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")
    def numberOfPermutations = 100
    def test_dir = new File("out/test").mkdir()
    Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")


    def "test uniform grid"()
    {
        given:
        def outputFileName = new File("out/test/test_uniform.csv")

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
        cmi closeTo(0.09727, 0.09728)
    }

    def "test pseudouniform grid"()
    {
        given:
        def outputFileName = new File("out/test/test_pseudouniform.csv")

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
        cmi closeTo(0.09727, 0.09728)
    }

    def "test cupid"()
    {
        given:
        def outputFileName = new File("out/test/test_cupid.csv")

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
        cmi closeTo(0.09727, 0.09728)
    }

    def "test normal"()
    {
        given:
        def outputFileName = new File("out/test/test_iterative_partitioning.csv")

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
        cmi closeTo(0.09727, 0.09728)
    }

    def "test normal with set format"()
    {
        given:
        genesMiRNA = new File("data/10_genes_mirna_interactions_set_format.txt")
        def outputFileName = new File("out/test/test_iterative_partitioning_set.csv")

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
        cmi closeTo(0.09727, 0.09728)
    }

}