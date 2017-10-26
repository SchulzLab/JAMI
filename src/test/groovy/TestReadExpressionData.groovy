import org.mpii.jami.ExpressionData
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestReadExpressionData extends Specification {

    def "read full gene expression table"(){
        given:
        def red = new ExpressionData()
        def fileGeneExpr = new File("data/10_genes_gene_expr.txt")

        when:
        red.readFile(fileGeneExpr, true)

        then:
        red.getExpressionData().size() == 10
        red.getExpressionData().get(0).size() == 362
    }

    def "read full mir expression table"(){
        given:
        def red = new ExpressionData()
        def fileGeneExpr = new File("data/10_genes_mir_expr.txt")

        when:
        red.readFile(fileGeneExpr, true)

        then:
        red.getExpressionData().size() == 79
        red.getExpressionData().get(0).size() == 362
    }

    def "read only three records from table"(){
        given:
        def red = new ExpressionData(["ENSG00000110427", "ENSG00000151746", "ENSG00000135631"])
        def fileGeneExpr = new File("data/10_genes_gene_expr.txt")

        when:
        red.readFile(fileGeneExpr, true)

        then:
        red.getExpressionData().size() == 3
        red.getExpressionData().get(0).size() == 362
    }

    def "read file without header"(){
        given:
        def red = new ExpressionData()
        def fileGeneExpr = new File("data/single_gene_pair_mir_expr.txt")

        when:
        red.readFile(fileGeneExpr, false)

        then:
        red.getExpressionData().size() == 6
        red.getExpressionData().get(0).size() == 362
    }
}