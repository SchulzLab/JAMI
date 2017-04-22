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
    private ArrayList<String[]> listOfTriples;
    private boolean tripleFormat;
    public static int numberOfSamples=364;  //this can be removed for the transposed format

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
        geneExpr=new ReadExpressionData(geneNames);



        listOfTriples=new ArrayList<>();

        for (int i = 0; i < geneNames.size(); i++) {
            String gene=geneNames.get(i);
            HashSet<String> gene1set=new HashSet<>();
            for (int j = 0; j < genesToMiRNA.get(gene).length; j++) {
                gene1set.add(genesToMiRNA.get(gene)[j]);
            }
            for (int j = 0; j < i; j++) {
                String gene2=geneNames.get(j);
                String[] gene2list=genesToMiRNA.get(gene2);
                for (String miRNAid : gene2list) {
                    if(gene1set.contains(miRNAid)){
                        String[] triple={gene,gene2,miRNAid};
                        listOfTriples.add(triple);
                    }
                }

            }
        }
    }

    private void triplesVariant2(){
       listOfTriples= ReadFiles.readTriples(fileGenesMiRNA, "\t", false);
        HashSet<String> geneNamesHS=new HashSet<>();
        HashSet<String> namesMiRNAHS=new HashSet<>();

        int counter=0;
        int bigCounter=0;
        for (String[] stringTriple : listOfTriples) {
            geneNamesHS.add(stringTriple[0]);
            geneNamesHS.add(stringTriple[1]);
            namesMiRNAHS.add(stringTriple[2]);
//            counter++;
//            bigCounter++;
//            if(counter==100){
//                counter=0;
//                System.out.println("big counter: "+bigCounter);
//            }
        }

        ArrayList<String> geneNames=new ArrayList<>(geneNamesHS);
        ArrayList<String> namesMiRNA=new ArrayList<>(namesMiRNAHS);

        miRExpr=new ReadExpressionData(namesMiRNA);
        geneExpr=new ReadExpressionData(geneNames);

       System.out.println("after file read");


        System.out.println("String to int");


    }

    public void runComputation(){

        System.out.println("start");
        if(tripleFormat){
            triplesVariant2();
        }
        else{
            triplesVariant1();
        }

        System.out.println("triples read");
      //  geneExpr.numberOfSamples =20531;

        geneExpr.numberOfSamples =this.numberOfSamples; //the three lines are to be commented for the new format
        geneExpr.genesAreRows=false;
        geneExpr.skipFirst=true;

      //  geneExpr.genesAreRows=true;
       // geneExpr.skipFirst=false;


        geneExpr.separator=this.separator;
        geneExpr.modifyName=false;
        geneExpr.readFile(fileGeneExpr);


        System.out.println("genes read");

        //miRExpr.separator=separator;

        miRExpr.numberOfSamples=this.numberOfSamples;   //the three lines are to be commented for the new format
        miRExpr.genesAreRows=false;
        miRExpr.skipFirst=true;

//        miRExpr.genesAreRows=true;
//        miRExpr.skipFirst=false;

        miRExpr.separator=this.separator;  //to be changed to \t in final
        miRExpr.modifyName=false;
        miRExpr.readFile(filemiRExpr);

        System.out.println("miRNA read");

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

            for (String[] oneTriple : listOfTriples) {
                Integer gene1Index=geneExpr.getNameToData().get(oneTriple[0]);
                Integer gene2Index=geneExpr.getNameToData().get(oneTriple[1]);
                Integer miRNAIndex=miRExpr.getNameToData().get(oneTriple[2]);

                if(gene1Index!=null&&gene2Index!=null&&miRNAIndex!=null){
                    computeCMIAndStore(gene1Index,gene2Index,miRNAIndex);
                }
                else{
                    System.out.println("not existing name");
                }

            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }



        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - timeStart));


    }

    //TODO work with the triples only by means of strings. Expression file structures are reduced after expression reading, the indices are different!

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
