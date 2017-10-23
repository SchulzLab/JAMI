import org.mpii.jami.CompleteRun
import spock.lang.Specification


/**
 * Created by mlist on 10/23/17.
 */
class TestVariousInputFormats extends Specification {
    def "input triplets and gene/miRNA expression files with samples in rows without parallel"() {
        given:
        def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
        def fileGeneExpr = new File("data/125_genes_gene_expr_to_test_in_cupid.txt")
        def filemiRNAExpr = new File("data/125_genes_mir_expr_to_test_in_cupid.txt")
        def outputFileName = new File("test_triplets_no_parallel_samples_in_rows.csv")
        def numberOfPermutations = 1000
        def separator = "\t"
        def parallel = false
        def numberOfSamples = 50
        def tripleFormat = true

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,parallel,separator,tripleFormat,numberOfSamples);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "input triplets and gene/miRNA expression files with samples in rows with parallel execution"() {
        given:
        def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
        def fileGeneExpr = new File("data/125_genes_gene_expr_to_test_in_cupid.txt")
        def filemiRNAExpr = new File("data/125_genes_mir_expr_to_test_in_cupid.txt")
        def outputFileName = new File("test_triplets_parallel_samples_in_rows.csv")
        def numberOfPermutations = 1000
        def separator = "\t"
        def parallel = true
        def numberOfSamples = 50
        def tripleFormat = true

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,parallel,separator,tripleFormat,numberOfSamples);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }
}