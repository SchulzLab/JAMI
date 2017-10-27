import org.mpii.jami.CompleteRun
import spock.lang.Specification

/**
 * Created by mlist on 10/26/17.
 */
class TestPValueCutoff extends Specification {

    def genesMiRNA = new File("data/10_genes_mirna_interactions_triplet_format.txt")
    def fileGeneExpr = new File("data/10_genes_gene_expr.txt")
    def filemiRNAExpr = new File("data/10_genes_mir_expr.txt")
    def numberOfPermutations = 100

    def test_dir = new File("out/test").mkdir()

    def "test p-value cutoff"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("out/test/test_p_value_cutoff.csv")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat,
                "", 0, -1, true, 0.05);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk < 100
    }
}