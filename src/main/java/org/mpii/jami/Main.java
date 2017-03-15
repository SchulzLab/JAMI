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

        //RunExamples.oneGenePair("BRCA1","CDH1","data/test_gene_expr.txt","data/test_mirna_expr.txt", "data/test_outputCMI.txt");

        String fileGenePairs=args[0];
        String genesMiRNA=args[1];
        String fileGegeExpr=args[2];
        String filemiRNAExpr=args[3];

//        String fileGenePairs="/home/fuksova/First Project/new_data/backup/genePairs.txt"; //I need to check on the final source files
//        String genesMiRNA="/home/fuksova/First Project/new_data/backup/genesMirna.txt";   //In parallel version, the p value is zero
//        String fileGegeExpr="/home/fuksova/First Project/new_data/backup/Gene_expr_new.txt";
//        String filemiRNAExpr="/home/fuksova/First Project/new_data/backup/miRNA_expr_liver";

        boolean parallel=true;
        String outputFileName="outputPar.csv";
        String separator="\t";
        int numberOfPermutations=1000;


        CompleteRun completeRun=new CompleteRun(fileGenePairs,genesMiRNA,fileGegeExpr,filemiRNAExpr,outputFileName,numberOfPermutations,separator,parallel);
        completeRun.runComputation();
    }









}


