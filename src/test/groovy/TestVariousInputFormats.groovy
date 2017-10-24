import org.mpii.jami.CompleteRun
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestVariousInputFormats extends Specification {

    def numberOfPermutations = 100
    def separator = "\t"

    def "input triplets and gene/miRNA expression files with samples in rows"() {
        given:
        def numberOfSamples = 50
        def tripleFormat = true
        def outputFileName = new File("test_triplets_samples_in_rows.csv")
        def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
        def fileGeneExpr = new File("data/125_genes_gene_expr_to_test_in_cupid.txt")
        def filemiRNAExpr = new File("data/125_genes_mir_expr_to_test_in_cupid.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,numberOfSamples);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "input triplets and gene/miRNA expression files with samples in columns"() {
        given:
        def tripleFormat = true
        def outputFileName = new File("test_triplets_samples_in_rows.csv")
        def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
        def fileGeneExpr = new File("data/125_genes_gene_expr_sample_in_cols.txt")
        def filemiRNAExpr = new File("data/125_genes_mir_expr_sample_in_cols.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,separator,tripleFormat,0);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

}