package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import org.mpii.jami.model.Triplet
import spock.lang.Shared
import spock.lang.Specification

import static spock.util.matcher.HamcrestMatchers.closeTo


/**
 * Created by mlist on 12/1/17.
 */
class TestDealingWithZeros extends Specification {

    @Shared
    def genesMiRNA = new File("data/single_gene_pair_triplets.txt")

    @Shared
    def fileGeneExpr = new File("data/single_gene_pair_gene_expr.txt")

    @Shared
    def filemiRNAExpr = new File("data/single_gene_pair_mir_expr.txt")

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

    def "test one gene pair with zeros"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_basics.txt")
        File modulatorGenewithZeros = new File("data/single_gene_pair_gene_expr_61_zeros.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,modulatorGenewithZeros,filemiRNAExpr, outputFileName);
        completeRun.setHeader(false);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 2
        Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double cmi = t.getCmi()
        cmi closeTo(0.1461176, 0.00001)

        def query_reverse = new Triplet("ENSG00000105855", "ENSG00000100767", "MIMAT0000421")
        def t_reverse = completeRun.getTriplets().find{
            (it == query_reverse)
        }
        def cmi_reverse = t_reverse.getCmi()
        cmi_reverse closeTo(0.0911841, 0.00001)
    }

    def "test one gene pair miRNA is zero"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_basics_mirna_expr_zero.txt")
        File miRNAwithZeros = new File("data/single_gene_pair_mir_expr_zeros.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,miRNAwithZeros, outputFileName);
        completeRun.setHeader(false);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 2
        Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double pValue = t.getpValue()
        pValue closeTo(1, 0.05)
    }

    def "test one gene pair modulator gene is zero"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_basics_gene_expr_zero.txt")
        File modulatorGenewithZeros = new File("data/single_gene_pair_gene_expr_zeros.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,modulatorGenewithZeros,filemiRNAExpr, outputFileName);
        completeRun.setHeader(false);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 2
        Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double pValue = t.getpValue()
        pValue closeTo(1, 0.05)

        Triplet reverseQuery = new Triplet("ENSG00000105855", "ENSG00000100767", "MIMAT0000421")
        Triplet reverseT = completeRun.getTriplets().find{
            (it == reverseQuery)
        }
        double reversePValue = reverseT.getpValue()
        reversePValue closeTo(1, 0.05)
    }


    def "test one gene pair modulator gene and miRNA is zero"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_basics_gene_and_mir_expr_zero.txt")
        File modulatorGenewithZeros = new File("data/single_gene_pair_gene_expr_zeros.txt")
        File miRNAwithZeros = new File("data/single_gene_pair_mir_expr_zeros.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,modulatorGenewithZeros,miRNAwithZeros, outputFileName);
        completeRun.setHeader(false);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 2
        Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")
        Triplet t = completeRun.getTriplets().find{
            (it == query)
        }
        double pValue = t.getpValue()
        pValue closeTo(1, 0.05)

        Triplet reverseQuery = new Triplet("ENSG00000105855", "ENSG00000100767", "MIMAT0000421")
        Triplet reverseT = completeRun.getTriplets().find{
            (it == reverseQuery)
        }
        double reversePValue = reverseT.getpValue()
        reversePValue closeTo(1, 0.05)
    }

}