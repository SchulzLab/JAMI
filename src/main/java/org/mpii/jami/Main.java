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
//
//        boolean tripleFormat=true;
//        String genesMiRNA="/home/fuksova/MMCI/data_3_17/125_genes_ceRNA_interactions_to_test_in_cupid.txt";
//        String filemiRNAExpr="/home/fuksova/MMCI/data_3_17/125_genes_mir_expr_to_test_in_cupid.txt";
//        String fileGeneExpr="/home/fuksova/MMCI/data_3_17/125_genes_gene_expr_to_test_in_cupid.txt";
//
//        String outputFileName="outputPar.csv";
//        String separator="\t";
//        int numberOfPermutations=1000;
//
//        int numberOfSamples=364;
//        boolean parallel=false;


        boolean tripleFormat=true;
        String genesMiRNA;
        String filemiRNAExpr;
        String fileGeneExpr;
        String outputFileName;
        int numberOfSamples;

        String separator="\t";
        int numberOfPermutations=1000;
        boolean parallel=false;


        if(args.length>=4&&args.length<=6){
            genesMiRNA=args[0];
            fileGeneExpr=args[1];
            filemiRNAExpr=args[2];
            outputFileName=args[3];
            if(args.length>=5){
                if(args[4].equals("t")){
                    tripleFormat=true;
                }
                else if(args[4].equals("s")){
                    tripleFormat=false;
                }
                else{
                    throw new IllegalArgumentException("Invalid settings for interaction file. Use either 's' or 't'.");
                }
            }
            else{
                tripleFormat=true;
            }
            if(args.length==6){
                numberOfSamples=Integer.valueOf(args[5]);
            }
            else{
                numberOfSamples=-1;
            }


            CompleteRun completeRun=new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,outputFileName,numberOfPermutations,parallel,separator,tripleFormat,numberOfSamples);
            completeRun.runComputation();
        }
        else{
            System.out.println("Input format:");
            System.out.println("First four arguments are compulsory: 0-File with interactions, 1-File with gene expression, 2-File with miRNA expression,3-Output file name.");
            System.out.println("4-Type of interaction file: 's' for set format, 't' for triple format, triple format is default.");
            System.out.println("5-Number of samples. If this is set, column format is assumed, otherwise row format is assumed.");
        }



    }









}


