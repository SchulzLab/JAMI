package org.mpii.jami.test

import org.mpii.jami.input.InteractionData
import org.mpii.jami.model.Triplet
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
        interactions.getTriplets().size() == 342 //number of triplets in file (or up to * 2 if reverse are not included)
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
        interactions.getTriplets().size() == 342
    }

    def "read interactions in gzipped set format"(){
        given:
        def interactions = new InteractionData()
        def tripletFile = new File("data/10_genes_mirna_interactions_set_format.txt.gz")

        when:
        interactions.readFileInSetFormat(tripletFile)

        then:
        interactions.getGenes().size() == 10 //8 genes have no interactions
        interactions.getMiRNAs().size() == 32
        interactions.getTriplets().size() == 342
    }

    def "read large interaction file in gzipped set format"(){
        given:
        def interactions = new InteractionData()
        def tripletFile = new File("data/mircode_set_format.txt.gz")

        when:
        interactions.readFileInSetFormat(tripletFile, "ENSG00000171862")

        then:
        interactions.getGenes().size() == 50186
        interactions.getMiRNAs().size() == 130
        interactions.getTriplets().contains(new Triplet("ENSG00000171862", "ENSG00000237984", "MIMAT0000076")) == true
    }
}