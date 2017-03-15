package org.mpii.jami;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 3/14/17.
 */
public class CompleteRun {
    private String outputFileName;
    private String fileGenePairs;
    private String fileGenesMiRNA;
    private  String fileGeneExpr;
    private String filemiRExpr;
    private boolean parallel;
    private int numberOfPermutations;
    private ReadExpressionData geneExpr;
    private ReadExpressionData miRExpr;
    private  BufferedWriter bw;
    private String separator;
    HashMap<String,ArrayList<Integer>> genesToMiRNAInt;

    public CompleteRun(String fileGenePairs,String fileGenesMiRNA,String fileGeneExpr,String filemiRExpr,String outputFileName,int numberOfPermutations,String separator,boolean parallel){
        this.outputFileName=outputFileName;
        this.fileGenePairs=fileGenePairs;
        this.fileGenesMiRNA=fileGenesMiRNA;
        this.fileGeneExpr=fileGeneExpr;
        this.filemiRExpr=filemiRExpr;
        this.parallel=parallel;
        this.separator=separator;

        this.numberOfPermutations=numberOfPermutations;
    }


    public void runComputation(){
        ArrayList<String[]> listPairs=new ArrayList<>();
        String pairSplit=separator;
        Pattern pattern = Pattern.compile(pairSplit);
        ArrayList<String> geneNames=new ArrayList<>();
        try {
            BufferedReader br=new BufferedReader(new FileReader(fileGenePairs));
            String line = br.readLine();
            while(line!=null){
                String[] pair = pattern.split(line);
                listPairs.add(pair);
                geneNames.add(pair[0]);
                geneNames.add(pair[1]);
                line=br.readLine();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


        geneExpr=new ReadExpressionData(geneNames);

        geneExpr.numberOfSamples =20531; //the two lines are to be commented for the new format
        geneExpr.genesAreRows=true;
        geneExpr.skipFirst=false;

      //  geneExpr.skipFirst=true;
        geneExpr.separator="\t";
        //geneExpr.genesAreRows=false;
        geneExpr.modifyName=false;
        geneExpr.readFile(fileGeneExpr);

        ArrayList<String> namesMiRNA= ReadFiles.extractNamesFromLine(filemiRExpr," ",0);
        miRExpr=new ReadExpressionData(namesMiRNA);
        miRExpr.skipFirst=true;
        //miRExpr.separator=separator;
        miRExpr.separator=" ";  //to be changed to \t in final
        miRExpr.genesAreRows=false;
        miRExpr.modifyName=false;
        miRExpr.readFile(filemiRExpr);

        HashMap<String, Integer> miNameToData = miRExpr.getNameToData();

        HashMap<String, String[]> genesToMiRNA = ReadFiles.geneToMiRNA(fileGenesMiRNA, "\t", ",", false);
        genesToMiRNAInt=new HashMap<>();
        for (Map.Entry<String, String[]> stringEntry : genesToMiRNA.entrySet()) {
            int length = stringEntry.getValue().length;
            ArrayList<Integer> miRNAIntegers=new ArrayList<>();
            for (int i=0;i<length;i++){
                miRNAIntegers.add(miNameToData.get(stringEntry.getValue()[i]));
            }
            genesToMiRNAInt.put(stringEntry.getKey(),miRNAIntegers);

        }

        int numberOfPermutations=1000;

        long timeStart=System.currentTimeMillis();

        FileWriter fw = null;
        try {
            fw = new FileWriter(outputFileName);
            bw = new BufferedWriter(fw);
            bw.write("geneInteracting");
            bw.write(separator);
            bw.write("geneConditional");
            bw.write(separator);
            bw.write("miRNA");
            bw.write(separator);
            bw.write("CMI");
            bw.write(separator);
            bw.write("pValue\n");
            for (int i = 0; i < listPairs.size(); i++) {
                String gene1=listPairs.get(i)[0];
                String gene2=listPairs.get(i)[1];
                HashSet<Integer> setOfMiRNA=new HashSet<>();
                setOfMiRNA.addAll(genesToMiRNAInt.get(gene1));
                ArrayList<Integer> intersection=new ArrayList<>();
                ArrayList<Integer> listToGene2=genesToMiRNAInt.get(gene2);
                for (Integer integer : listToGene2) {
                    if(setOfMiRNA.contains(integer)){
                        intersection.add(integer);
                    }
                }
                computeCMIAndStore(geneExpr.getNameToData().get(gene1),geneExpr.getNameToData().get(gene2));
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }



        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - timeStart));


    }





    public void computeCMIAndStore(int gene1,
                                          int gene2){
        double[] gene1Data = geneExpr.getExpressionData().get(gene1);
        double[] gene2Data = geneExpr.getExpressionData().get(gene2);
        ArrayList<double[]> exprDatamiRna = miRExpr.getExpressionData();
        String gene1Name=geneExpr.getIntegersToNames().get(gene1);
        String gene2Name=geneExpr.getIntegersToNames().get(gene2);
        //File outputFile=new File(outputFileName);

            HashSet<Integer> setOfMiRNA=new HashSet<>();
            setOfMiRNA.addAll(genesToMiRNAInt.get(gene1Name));
            ArrayList<Integer> intersection=new ArrayList<>();
            ArrayList<Integer> listToGene2=genesToMiRNAInt.get(gene2Name);
            for (Integer integer : listToGene2) {
                if(setOfMiRNA.contains(integer)){
                    intersection.add(integer);
                }
            }

        if(!intersection.isEmpty()) {
            try {

                //CMI(gene1Name,miRNA|gene2Name)
                for (int i = 0; i < intersection.size(); i++) {
                    int index=intersection.get(i);
                    String miRNAName=miRExpr.getIntegersToNames().get(index);
                    bw.write(gene1Name+separator+gene2Name+separator+miRNAName+separator);
                    double[] result= computeTriple(gene2Data, gene1Data,  index);
                    bw.write(result[0]+separator+result[1]+"\n");
                    bw.flush();

                }

                //CMI(gene2Name,miRNA|gene1Name)
                for (int i = 0; i < intersection.size(); i++) {
                    int index=intersection.get(i);
                    String miRNAName=miRExpr.getIntegersToNames().get(index);
                    bw.write(gene2Name+separator+gene1Name+separator+miRNAName+separator);
                    double[] result= computeTriple(gene1Data, gene2Data,  index);
                    bw.write(result[0]+separator+result[1]+"\n");
                    bw.flush();

                }
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

    }

    private double[] computeTriple(double[] geneCondData, double[] geneInterData, int i) {  //gene1Data is randomized
        ArrayList<double[]> data = new ArrayList<>();
        double[] result=new double[2];
        data.add(geneInterData);   //last added is randomized
        data.add(miRExpr.getExpressionData().get(i));
        data.add(geneCondData);
        CMIComplete cmiComplete;
        //cmiComplete.cmiAndPValuePseudoUniform(8);
        if(!parallel){
            CMIComplete.initRandomized(geneInterData.length,numberOfPermutations);
            cmiComplete = new CMIComplete(numberOfPermutations, data);
            cmiComplete.computeCMIandPValueBetterPartitioning();
        }
        else{
            cmiComplete = new CMIComplete(numberOfPermutations, data, false);
            cmiComplete.computeCMIandPValueParallel();
        }
        result[0]=cmiComplete.cmi;
        result[1]=cmiComplete.pValue;
        return result;

    }




}
