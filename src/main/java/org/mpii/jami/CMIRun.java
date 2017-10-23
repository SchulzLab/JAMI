package org.mpii.jami;

import java.io.*;
import java.util.ArrayList;
import java.util.regex.Pattern;

/**
 * Created by mlist on 12/15/16.
 */
@SuppressWarnings("Duplicates")
public class CMIRun {

    public static void computeCMIAndStore(ReadExpressionData geneData,
                                          ReadExpressionData miRNAData,
                                          int gene1,
                                          int gene2,
                                          String outputFileName,
                                          int numberOfPermutations,
                                          boolean parallel){
        double[] gene1Data = geneData.getExpressionData().get(gene1);
        double[] gene2Data = geneData.getExpressionData().get(gene2);
        ArrayList<double[]> exprDatamiRna = miRNAData.getExpressionData();
        String gene1Name=geneData.getIntegersToNames().get(gene1);
        String gene2Name=geneData.getIntegersToNames().get(gene2);
        //File outputFile=new File(outputFileName);
        FileWriter fw = null;
        try {
            fw = new FileWriter(outputFileName);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("Interaction of "+gene1Name+" with "+gene2Name+"\n");
            bw.write("I("+gene1Name+",miRNA|"+gene2Name+")\n");
            bw.write("miRNA id\tCMI\tp-value\n");
            for (int i = 0; i < exprDatamiRna.size(); i++) {
                writeOutput(miRNAData, numberOfPermutations, parallel, gene2Data, gene1Data, exprDatamiRna, bw, i);
            }
            bw.write("\nI("+gene2Name+",miRNA|"+gene1Name+")\n");
            bw.write("miRNA id\tCMI\tp-value\n");
            for (int i = 0; i < exprDatamiRna.size(); i++) {
                writeOutput(miRNAData, numberOfPermutations, parallel, gene1Data, gene2Data, exprDatamiRna, bw, i);
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void writeOutput(ReadExpressionData miRNAData, int numberOfPermutations, boolean parallel, double[] gene1Data, double[] gene2Data, ArrayList<double[]> exprDatamiRna, BufferedWriter bw, int i) throws IOException {
        ArrayList<double[]> data = new ArrayList<>();
        data.add(gene2Data);
        data.add(exprDatamiRna.get(i));
        data.add(gene1Data);
        CMIComplete cmiComplete;
        //cmiComplete.cmiAndPValuePseudoUniform(8);
        if(!parallel){
            cmiComplete = new CMIComplete(numberOfPermutations, data);
            cmiComplete.computeCMIandPValueBetterPartitioning();
        }
        else{
            cmiComplete = new CMIComplete(numberOfPermutations, data, false);
            cmiComplete.computeCMIandPValueParallel();
        }
        bw.write(miRNAData.getIntegersToNames().get(i)+"\t"+cmiComplete.cmi+"\t"+cmiComplete.pValue+"\n");
        bw.flush();
    }

}
