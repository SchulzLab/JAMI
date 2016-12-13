package org.mpii.jami;

/**
 * Created by fuksova on 10/7/15.
 */
public class Main {

    /**
     * Uncomment one of the options
     * @param args
     */
    public static void main(String []args) {
        //CMI computations
        //RunExamples.readExpressionDataTest();
        //RunExamples.basicCupidExampleDemo();
        //RunExamples.oneGenePair("ABHD1","CDKN2B","new_data/Gene_expr_new.txt","new_data/miRNA_expr_liver","data/testOutputCMI.txt");
        //RunExamples.moreGenePairs("data/genePairs.txt", "new_data/Gene_expr_new.txt", "new_data/miRNA_expr_liver");

        //Score exraction and gene pairs selection
        //GenePairsSelection.inspectContingencyTables();
        //GenePairsSelection.prunnedPValuesComplete();
        //GenePairsSelection.qValues();
        //GenePairsSelection.rewritePValues();
        //GenePairsSelection.writeGenePairsBellowP();

        RunExamples.oneGenePair("BRCA1","CDH1","data/test_gene_expr.txt","data/test_mirna_expr.txt", "data/test_outputCMI.txt");
    }









}


