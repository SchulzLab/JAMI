import org.mpii.jami.CompleteRun
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestVariousCmiImplementations extends Specification {

    def genesMiRNA = new File("data/125_genes_ceRNA_interactions_to_test_in_cupid.txt")
    def fileGeneExpr = new File("data/125_genes_gene_expr_sample_in_cols.txt")
    def filemiRNAExpr = new File("data/125_genes_mir_expr_sample_in_cols.txt")
    def numberOfPermutations = 100

    def "test uniform grid"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("test_uniform.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat,
                "uniform", 3, -1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "test pseudouniform grid"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("test_pseudouniform.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat,
                "pseudouniform", 3, -1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "test cupid"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("test_cupid.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat,
                "cupid", 0, -1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

    def "test normal"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("test_iterative_partitioning.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat,
                "", 0, -1);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 3866
    }

}