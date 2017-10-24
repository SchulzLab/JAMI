import org.mpii.jami.ExpressionData
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestReadExpressionData extends Specification {

    def "read full table"(){
        given:
        def red = new ExpressionData()
        def fileGeneExpr = new File("data/125_genes_gene_expr_sample_in_cols.txt")

        when:
        red.readFile(fileGeneExpr)

        then:
        red.getExpressionData().size() == 125
    }

    def "read only three records from table"(){
        given:
        def red = new ExpressionData(["TIMD4", "GRIP2", "ILK"])
        def fileGeneExpr = new File("data/125_genes_gene_expr_sample_in_cols.txt")

        when:
        red.readFile(fileGeneExpr)

        then:
        red.getExpressionData().size() == 3
    }
}