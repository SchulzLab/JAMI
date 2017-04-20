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


        String genesMiRNA=args[0];
        String fileGegeExpr=args[1];
        String filemiRNAExpr=args[2];

//        String fileGenePairs="/home/fuksova/First Project/new_data/backup/genePairs.txt"; //I need to check on the final source files
//        String genesMiRNA="/home/fuksova/First Project/new_data/backup/genesMirna.txt";   //In parallel version, the p value is zero
//        String fileGegeExpr="/home/fuksova/First Project/new_data/backup/Gene_expr_new.txt";
//        String filemiRNAExpr="/home/fuksova/First Project/new_data/backup/miRNA_expr_liver";

        boolean parallel=false;
        String outputFileName="outputPar.csv";
        String separator="\t";
        int numberOfPermutations=1000;
        boolean tripleFormat=false;

        if(args.length>=4&&args[3].equals("t")){
            tripleFormat=true;
        }


        CompleteRun completeRun=new CompleteRun(genesMiRNA,fileGegeExpr,filemiRNAExpr,outputFileName,numberOfPermutations,separator,parallel,tripleFormat);
        completeRun.runComputation();
    }









}


