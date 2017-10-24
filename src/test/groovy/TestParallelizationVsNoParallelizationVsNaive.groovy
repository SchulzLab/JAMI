import org.mpii.jami.CompleteRun
import spock.lang.Specification


/**
 * Created by mlist on 10/23/17.
 */
class TestParallelizationVsNoParallelizationVsNaive extends Specification {


    def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
    def fileGeneExpr = new File("data/125_genes_gene_expr_to_test_in_cupid.txt")
    def filemiRNAExpr = new File("data/125_genes_mir_expr_to_test_in_cupid.txt")
    def numberOfPermutations = 1000
    def separator = "\t"

    def "iterative partitioning with 1 core"() {
        given:
        def numberOfSamples = 50
        def tripleFormat = true
        def outputFileName = new File("test_triplets_no_parallel_samples_in_rows.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,numberOfSamples,
                "", 0, 1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "iterative partitioning with 2 cores"() {
        given:
        def numberOfSamples = 50
        def tripleFormat = true
        def outputFileName = new File("test_triplets_parallel_samples_in_rows.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,numberOfSamples,
                "", 0, 2);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "iterative Partitioning -1 cores"() {
        given:
        def numberOfSamples = 50
        def tripleFormat = true
        def outputFileName = new File("test_triplets_parallel_samples_in_rows.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,numberOfSamples);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "1 core and CUPID implementation"() {
        given:
        def numberOfSamples = 50
        def tripleFormat = true
        def outputFileName = new File("test_triplets_parallel_samples_in_rows.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,numberOfSamples,
                "cupid",0,1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

}