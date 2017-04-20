package org.mpii.jami;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 3/14/17.
 */
public class CompleteRun {
    private String outputFileName;
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
    private ArrayList<int[]> listOfTriples;
    private boolean tripleFormat;

    public CompleteRun(String fileGenesMiRNA,String fileGeneExpr,String filemiRExpr,String outputFileName,int numberOfPermutations,String separator,boolean parallel,boolean tripleFormat){
        this.outputFileName=outputFileName;
        this.fileGenesMiRNA=fileGenesMiRNA;
        this.fileGeneExpr=fileGeneExpr;
        this.filemiRExpr=filemiRExpr;
        this.parallel=parallel;
        this.separator=separator;
        this.tripleFormat=tripleFormat;
        this.numberOfPermutations=numberOfPermutations;
    }


    private void triplesVariant1(){
        HashMap<String, String[]> genesToMiRNA = ReadFiles.geneToMiRNA(fileGenesMiRNA, "\t", ",", false);

        ArrayList<String> geneNames=new ArrayList<>();
        geneNames.addAll(genesToMiRNA.keySet());
        HashSet<String> miRNAnamesToUse=new HashSet<>();

        for (String[] miRNANames : genesToMiRNA.values()) {
            for (int i = 0; i < miRNANames.length; i++) {
                miRNAnamesToUse.add(miRNANames[i]);
            }
        }

        ArrayList<String> namesMiRNA=new ArrayList<>(miRNAnamesToUse);
        miRNAnamesToUse.clear();
        miRExpr=new ReadExpressionData(namesMiRNA);
        HashMap<String, Integer> miNameToData = miRExpr.getNameToData();
        geneExpr=new ReadExpressionData(geneNames);

        genesToMiRNAInt=new HashMap<>();
        for (Map.Entry<String, String[]> stringEntry : genesToMiRNA.entrySet()) {
            int length = stringEntry.getValue().length;
            ArrayList<Integer> miRNAIntegers=new ArrayList<>();
            for (int i=0;i<length;i++){
                miRNAIntegers.add(miNameToData.get(stringEntry.getValue()[i]));
            }
            genesToMiRNAInt.put(stringEntry.getKey(),miRNAIntegers);
        }

        listOfTriples=new ArrayList<>();

        for (int i = 0; i < geneNames.size(); i++) {
            String gene=geneNames.get(i);
            int gene1Index=geneExpr.getNameToData().get(gene);
            HashSet<Integer> gene1set=new HashSet<>(genesToMiRNAInt.get(gene));
            for (int j = 0; j < i; j++) {
                String gene2=geneNames.get(j);
                int gene2index=geneExpr.getNameToData().get(gene2);
                ArrayList<Integer> gene2list=genesToMiRNAInt.get(gene2);
                for (Integer miRNAid : gene2list) {
                    if(gene1set.contains(miRNAid)){

                        int[] triple={gene1Index,gene2index,miRNAid};
                        listOfTriples.add(triple);
                    }
                }

            }
        }
    }

    private void triplesVariant2(){
        ArrayList<String[]> stringTriples = ReadFiles.readTriples(fileGenesMiRNA, "\t", false);
        HashSet<String> geneNamesHS=new HashSet<>();
        HashSet<String> namesMiRNAHS=new HashSet<>();
        for (String[] stringTriple : stringTriples) {
            geneNamesHS.add(stringTriple[0]);
            geneNamesHS.add(stringTriple[1]);
            namesMiRNAHS.add(stringTriple[2]);
        }

        ArrayList<String> geneNames=new ArrayList<>(geneNamesHS);
        ArrayList<String> namesMiRNA=new ArrayList<>(namesMiRNAHS);

        for (String[] stringTriple : stringTriples) {
            int index1=geneExpr.getNameToData().get(stringTriple[0]);
            int index2=geneExpr.getNameToData().get(stringTriple[1]);
            int index3=miRExpr.getNameToData().get(stringTriple[2]);
            int[] oneTriple={index1,index2,index3};
            listOfTriples.add(oneTriple);
        }

        miRExpr=new ReadExpressionData(namesMiRNA);
        geneExpr=new ReadExpressionData(geneNames);
    }

    public void runComputation(){

        if(tripleFormat){
            triplesVariant2();
        }
        else{
            triplesVariant1();
        }

        geneExpr.numberOfSamples =20531; //the three lines are to be commented for the new format
        geneExpr.genesAreRows=true;
        geneExpr.skipFirst=false;

      //  geneExpr.skipFirst=true;
        geneExpr.separator="\t";
        //geneExpr.genesAreRows=false;
        geneExpr.modifyName=false;
        geneExpr.readFile(fileGeneExpr);


        miRExpr.skipFirst=true;
        //miRExpr.separator=separator;
        miRExpr.separator=" ";  //to be changed to \t in final
        miRExpr.genesAreRows=false;
        miRExpr.modifyName=false;
        miRExpr.readFile(filemiRExpr);


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

            for (int[] oneTriple : listOfTriples) {
                computeCMIAndStore(oneTriple[0],oneTriple[1],oneTriple[2]);
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }



        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - timeStart));


    }


    public void computeCMIAndStore(int gene1, int gene2, int miRNA) {
        double[] gene1Data = geneExpr.getExpressionData().get(gene1);
        double[] gene2Data = geneExpr.getExpressionData().get(gene2);
        double[] miRNAData=miRExpr.getExpressionData().get(miRNA);
        String gene1Name = geneExpr.getIntegersToNames().get(gene1);
        String gene2Name = geneExpr.getIntegersToNames().get(gene2);
        String miRNAName = miRExpr.getIntegersToNames().get(miRNA);


        try {

            bw.write(gene1Name + separator + gene2Name + separator + miRNAName + separator);
            double[] result = computeTriple(gene2Data, gene1Data, miRNAData);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();

            bw.write(gene2Name + separator + gene1Name + separator + miRNAName + separator);
            result = computeTriple(gene1Data, gene2Data, miRNAData);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();


        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private double[] computeTriple(double[] geneCondData, double[] geneInterData, double[] miRNAData) {  //gene1Data is randomized
        ArrayList<double[]> data = new ArrayList<>();
        double[] result=new double[2];
        data.add(geneInterData);   //last added is randomized
        data.add(miRNAData);
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
