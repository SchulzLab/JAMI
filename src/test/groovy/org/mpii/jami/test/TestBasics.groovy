package org.mpii.jami.test

import org.mpii.jami.CompleteRun
import org.mpii.jami.cmi.Cube
import org.mpii.jami.helpers.AdditionalComputations
import org.mpii.jami.helpers.BenjaminiHochberg
import org.mpii.jami.helpers.FisherMethod
import org.mpii.jami.model.Triplet
import spock.lang.Shared
import spock.lang.Specification
import static spock.util.matcher.HamcrestMatchers.closeTo

/**
 * Created by mlist on 10/26/17.
 */
class TestBasics extends Specification {

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

    def "test one gene pair"()
    {
        given:
        def outputFileName = new File(testDir.absolutePath + "/test_basics.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr, outputFileName);
        //completeRun.setMethod("cupid")
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
        cmi closeTo(0.098, 0.00001)
    }

    def "test adjusting p-values"()
    {
        given:
        double[] pValues = [0.01, 0.0004, 0.000005, 0.0004, 0.5, 0.1, 0.000001, 1.0]

        when:
        double[] adjustedPValues = BenjaminiHochberg.adjustPValues(pValues)
        adjustedPValues = adjustedPValues.collect{ it -> Math.round(it * 10000000d) / 10000000d }

        then:
        adjustedPValues == [0.0160000, 0.0008000, 0.0000200, 0.0008000, 0.5714286, 0.1333333, 0.0000080, 1.0000000]
    }


    def "test Fisher method"()
    {
        given:
        List<Double> pValues =  [0.0001d, 0.0001d, 0.9999d, 0.9999d]

        when:
        double metaP = FisherMethod.combinedPValue(pValues)

        then:
        metaP closeTo(0.0000123, 0.0000001)
    }
}