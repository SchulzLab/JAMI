import org.mpii.jami.ExpressionData
import org.mpii.jami.InteractionData
import spock.lang.Specification


/**
 * Created by mlist on 10/24/17.
 */
class TestReadInteractionsData extends Specification {
    def "read interactions in triple format"(){
        given:
        def interactions = new InteractionData()
        def tripletFile = new File("data/10_genes_mirna_interactions_triplet_format.txt")

        when:
        interactions.readFileWithTriples(tripletFile)

        then:
        interactions.getGenes().size() == 10
        interactions.getMiRNAs().size() == 32 //number of miRNAs
        interactions.getTriplets().size() == 171 //number of triplets in file
    }

    def "read interactions in set format"(){
        given:
        def interactions = new InteractionData()
        def tripletFile = new File("data/10_genes_mirna_interactions_set_format.txt")

        when:
        interactions.readFileInSetFormat(tripletFile)

        then:
        interactions.getGenes().size() == 10 //8 genes have no interactions
        interactions.getMiRNAs().size() == 32
        interactions.getTriplets().size() == 171
    }
}