import org.mpii.jami.CompleteRun
import org.mpii.jami.model.Triplet
import spock.lang.Specification
import static spock.util.matcher.HamcrestMatchers.closeTo

/**
 * Created by mlist on 10/26/17.
 */
class TestBasics extends Specification {

    def genesMiRNA = new File("data/single_gene_pair_triplets.txt")
    def fileGeneExpr = new File("data/single_gene_pair_gene_expr.txt")
    def filemiRNAExpr = new File("data/single_gene_pair_mir_expr.txt")
    def numberOfPermutations = 1000

    def "test one gene pair"()
    {
        given:
        def tripleFormat = true
        def outputFileName = new File("out/test/test_basics.txt")

        when:
        CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,
                outputFileName,numberOfPermutations,tripleFormat, false);
        completeRun.runComputation();

        then:
        completeRun.completed == true
        completeRun.tripletsWrittenToDisk == 2
        Triplet query = new Triplet("ENSG00000100767", "ENSG00000105855", "MIMAT0000421")
        double cmi = (double) completeRun.getCmis().get(query)
        cmi closeTo(0.09727, 0.09728)
    }


}