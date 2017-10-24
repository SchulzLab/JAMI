import org.mpii.jami.ExpressionData
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestReadInteractionsData extends Specification {
    def "read full table"(){
        given:
        def red = new ExpressionData()
        def fileGeneExpr = new File("data/125_genes_gene_expr_sample_in_cols.txt")

        when:
        red.readFile(fileGeneExpr)

        then:
        red.getExpressionData().size() == 125
    }
}