package org.mpii.jami;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 1/19/16.
 * Some "demo" examples for demonstrating how to use the functions for CMI computation, expression data read etc.
 */
public class RunExamples {


    /**
     * This function demonstrates how to read files with gene expression data and how to access the obtained values
     * It just reads the data. It does not have any output
     */
    public static void readExpressionDataTest(){
        //Current version assumes that names of genes that we want to explore are given
        ArrayList<String> allNames=ReadFiles.extractNamesFromColumn("new_data/uniqueGeneList.txt", "\t", false);

        ReadExpressionData readColumn=new ReadExpressionData(allNames);
        readColumn.modifyName=true;
        readColumn.separator="\t";
        readColumn.skipFirst=true;
        readColumn.multiplier=1000000;

        //Reads the file with expression data
        readColumn.readFile("new_data/Gene_expr_new.txt");


        //Every gene name corresponds to certain integer which denotes index in expressionData where the expression
        //data for this gene can be found
        HashMap<String, Integer> namesToData = readColumn.getNameToData();
        //gets all expression data
        ArrayList<double[]> expressionData = readColumn.getExpressionData();
        double[] esr1Expression = readColumn.getExpressionDataToID("ESR1");

        ReadExpressionData readRow=new ReadExpressionData(allNames);
        readRow.modifyName=true;
        readRow.separator=" ";
        //readRow.numberOfSamples=20531;
        readRow.skipFirst=true;
        readRow.genesAreRows=false;
        readRow.multiplier=1000000;
        readRow.readFile("new_data/Genes_expr_liver");

    }


    /**
     * This method enables to run the basic Cupid example from cupidcerna.m. If certain rows are uncommented,
     * it is possible to use other CMI computation methods and compare p-values obtained for all triplets computed by
     * those various methods.
     */
    public static void newScoreComputationBetterIP(){
        ArrayList<String> allNames=ReadFiles.extractNamesFromColumn("data/expr1", "\t", false);
        ArrayList<String> geneNames=new ArrayList<>(2);
        geneNames.add(allNames.remove(0));
        geneNames.add(allNames.remove(0));
        Collections.sort(allNames);
        ReadExpressionData readGeneData=new ReadExpressionData(geneNames);
        readGeneData.separator="\t";
        readGeneData.skipFirst=false;
        readGeneData.modifyName=false;
        readGeneData.readFile("data/expr1");
         ArrayList<double[]> genesData=readGeneData.getExpressionData();


        ScoreStructure.geneInd=1;
        ScoreStructure.miRNAInd=2;
        ScoreStructure.valueInd=0;
        HashMap<String, HashMap<String, Double>> hmESR1 = ScoreStructure.readPosAndNegData("data/ESR1Pos.txt", "data/ESR1Neg.txt",false,false);
        HashMap<String, HashMap<String, Double>> hmHIF1A = ScoreStructure.readPosAndNegData("data/HIF1APos.txt", "data/HIF1ANeg.txt",false,false);

        ScoreTable st=new ScoreTable(hmESR1.get("ESR1"),hmHIF1A.get("HIF1A"),allNames,0,allNames.size());
        double cutoff=0.05;


        SelectMiRNACandidates pt=new SelectMiRNACandidates(st.getMiRNAnames(),1000,st.data);
        pt.selectCandidates(cutoff);
        int[] candidateIndices = pt.getSelectedIndices();


        ReadExpressionData readmiRNAData=new ReadExpressionData(allNames);
        readmiRNAData.separator="\t";
        readmiRNAData.skipFirst=false;
        readmiRNAData.modifyName=false;
        readmiRNAData.readFile("data/expr1");
        ArrayList<double[]> miRNAData=readmiRNAData.getExpressionData();

        double[] pvOrig1=new double[pt.getNumberSelected()];
        double[] pvChanged1=new double[pt.getNumberSelected()];

        double[] pvOrig2=new double[pt.getNumberSelected()];
        double[] pvChanged2=new double[pt.getNumberSelected()];


        long start=System.currentTimeMillis();
        for (int i = 0; i < pt.getNumberSelected(); i++) {
            //In file expr1 are stored expression data for two genes and for a set of miRNA. The data for genes
            //are at the begininning of the file
            ArrayList<double[]> data=new ArrayList<>();
            data.add(genesData.get(0));
            data.add(miRNAData.get(candidateIndices[i]));
            data.add(genesData.get(1));
            CMIComplete cmiComplete=new CMIComplete(1000,data);

            cmiComplete.cmiAndPValueFixedGrid();
            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue + ", ");
            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
            System.out.println(mediatorStatus);
            pvChanged1[i]=cmiComplete.pValue;

//            cmiComplete.computeCMIandPVUniformGrid(3);
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged1[i]=cmiComplete.cmi;
            //pvOrig1[i]=cmiComplete.pValue;

//            cmiComplete.cmiAndPValuePseudoUniform(4);
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged1[i]=cmiComplete.pValue;

//            cmiComplete.setMaxDeep(2);
//            cmiComplete.computeCMIandPValueBetterPartitioning();
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged1[i]=cmiComplete.pValue;

//            cmiComplete.setMaxDeep(Integer.MAX_VALUE);
//            cmiComplete.computeCMIandPValueBetterPartitioning();
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//           boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvOrig1[i]=cmiComplete.cmi;

        }
        for (int i = 0; i < pt.getNumberSelected(); i++) {
            ArrayList<double[]> data=new ArrayList<>();
            //otevirani souboru v kazdem behu je asi nevhyhodne a hlavne zpomaleni oproti orig. kodu. Nutno zmenit!
            data.add(genesData.get(1));                                         //V ramci priblizeni se orig. implementaci by to mozna chtelo nechat puvodni tabulku dat a pamatovat si indexy vybranych
            data.add(miRNAData.get(candidateIndices[i]));
            data.add(genesData.get(0));
            CMIComplete cmiComplete=new CMIComplete(1000,data);

            cmiComplete.cmiAndPValueFixedGrid();
            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
            System.out.println(mediatorStatus);
            pvChanged2[i]=cmiComplete.pValue;

//            cmiComplete.computeCMIandPVUniformGrid(3);
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged2[i]=cmiComplete.cmi;
            //pvOrig2[i]=cmiComplete.pValue;

//            cmiComplete.cmiAndPValuePseudoUniform(4);
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged2[i]=cmiComplete.pValue;

//            cmiComplete.setMaxDeep(2);
//            cmiComplete.computeCMIandPValueBetterPartitioning();
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvChanged2[i]=cmiComplete.pValue;

//            cmiComplete.setMaxDeep(Integer.MAX_VALUE);
//            cmiComplete.computeCMIandPValueBetterPartitioning();
//            System.out.print(pt.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue+", ");
//            boolean mediatorStatus=cmiComplete.mediatorStatus(pt.getSelectedWeights()[i],cutoff);
//            System.out.println(mediatorStatus);
//            pvOrig2[i]=cmiComplete.cmi;
        }
//        CompareDoubleArrays evaluate1=new CompareDoubleArrays(pvOrig1,pvChanged1);  //compare p-values obtained by two different methods
//        evaluate1.computeAndOutput();
//        CompareDoubleArrays evaluate2=new CompareDoubleArrays(pvOrig2,pvChanged2);
//        evaluate2.computeAndOutput();

        long end=System.currentTimeMillis();
        System.out.println("time: "+(end-start));

    }



    /**
     * The original example cupidcerna.m is replicated here but with more efficient implementation
     */
    public static void basicCupidExampleDemo() {
        //In file expr1 are stored expression data for two genes and for a set of miRNA. The data for genes
        //are at the begininning of the file
        ArrayList<String> allNames = ReadFiles.extractNamesFromColumn("data/expr1", "\t", false);
        ArrayList<String> geneNames = new ArrayList<>(2);
        geneNames.add(allNames.remove(0));
        geneNames.add(allNames.remove(0));
        Collections.sort(allNames);
        ReadExpressionData readGeneData = new ReadExpressionData(geneNames);
        readGeneData.separator="\t";
        readGeneData.skipFirst=false;
        readGeneData.modifyName=false;
        readGeneData.readFile("data/expr1");
        ArrayList<double[]> genesData = readGeneData.getExpressionData();


        ScoreStructure.geneInd = 1;
        ScoreStructure.miRNAInd = 2;
        ScoreStructure.valueInd = 0;
        HashMap<String, HashMap<String, Double>> hmESR1 = ScoreStructure.readPosAndNegData("data/ESR1Pos.txt", "data/ESR1Neg.txt",false,false);
        HashMap<String, HashMap<String, Double>> hmHIF1A = ScoreStructure.readPosAndNegData("data/HIF1APos.txt", "data/HIF1ANeg.txt",false,false);

        ScoreTable st = new ScoreTable(hmESR1.get("ESR1"), hmHIF1A.get("HIF1A"), allNames, 0, allNames.size());
        double cutoff = 0.05;


        SelectMiRNACandidates selectMiRNACandidates = new SelectMiRNACandidates(st.getMiRNAnames(), 1000, st.data);
        selectMiRNACandidates.selectCandidates(cutoff);
        int[] candidateIndices = selectMiRNACandidates.getSelectedIndices();


        ReadExpressionData readmiRNAData = new ReadExpressionData(allNames);
        readmiRNAData.separator="\t";
        readmiRNAData.skipFirst=false;
        readmiRNAData.modifyName=false;
        readmiRNAData.readFile("data/expr1");

        ArrayList<double[]> miRNAData = readmiRNAData.getExpressionData();


        long start = System.currentTimeMillis();
        for (int i = 0; i < selectMiRNACandidates.getNumberSelected(); i++) {
            ArrayList<double[]> data = new ArrayList<>();

            data.add(genesData.get(0));
            data.add(miRNAData.get(candidateIndices[i]));
            data.add(genesData.get(1));
            CMIComplete cmiComplete = new CMIComplete(1000, data);

            cmiComplete.computeCMIandPValueBetterPartitioning();
            System.out.print(selectMiRNACandidates.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue + ", ");

            boolean mediatorStatus = cmiComplete.mediatorStatus(selectMiRNACandidates.getSelectedWeights()[i], cutoff);
            System.out.println(mediatorStatus);
        }
        for (int i = 0; i < selectMiRNACandidates.getNumberSelected(); i++) {
            ArrayList<double[]> data = new ArrayList<>();
            data.add(genesData.get(1));
            data.add(miRNAData.get(candidateIndices[i]));
            data.add(genesData.get(0));
            CMIComplete cmiComplete = new CMIComplete(1000, data);


            cmiComplete.computeCMIandPValueBetterPartitioning();
            System.out.print(selectMiRNACandidates.getNames().get(candidateIndices[i]) + ": cmi: " + cmiComplete.cmi + ", p-value: " + cmiComplete.pValue + ", ");

            boolean mediatorStatus = cmiComplete.mediatorStatus(selectMiRNACandidates.getSelectedWeights()[i], cutoff);
            System.out.println(mediatorStatus);
        }
        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - start));
    }


    /**
     *
     * Computes complete CMI values for given pair of genes and a set of miRNA. Does not filter miRNA data.
     * @param gene1 name of first gene
     * @param gene2 name of second gene
     * @param fileGeneExpr File with gene expression data
     * @param filemiRExpr File with miRNA expression data
     * @param outputFileName Name of output file
     */
    public static void oneGenePair(String gene1,String gene2,String fileGeneExpr,String filemiRExpr,String outputFileName){
        ArrayList<String> geneNames=new ArrayList<>();
        geneNames.add(gene1);
        geneNames.add(gene2);
        ReadExpressionData geneExpr=new ReadExpressionData(geneNames);
        geneExpr.skipFirst=true;
        geneExpr.numberOfSamples=20531;
        geneExpr.separator="\t";
        geneExpr.genesAreRows=true;
        geneExpr.modifyName=true;
        geneExpr.readFile(fileGeneExpr);

        ArrayList<String> namesMiRNA= ReadFiles.extractNamesFromColumn(filemiRExpr,"\t", true);
        ReadExpressionData miRExpr=new ReadExpressionData(namesMiRNA);
        miRExpr.skipFirst=true;
        miRExpr.separator="\t";
        miRExpr.genesAreRows=true;
        miRExpr.modifyName=true;   //The names are extracted from the file including quotes, so the quotes should not be removed when comparing names
        miRExpr.readFile(filemiRExpr);

        int numberOfPermutations=1000;

        long timeStart=System.currentTimeMillis();
        computeCMIAndStore(geneExpr,miRExpr,gene1,gene2,outputFileName,numberOfPermutations);
        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - timeStart));
    }


    /**
     * Given a file with gene pairs and files with expression data, it computes CMI for every gene pair with all
     * miRNA and stores it into output files.
     * @param fileGenePairs File where gene pairs are stored in format: name of first gene \t name of second gene \n
     * @param fileGeneExpr File with gene expression data
     * @param filemiRExpr File with miRNA expression data
     */
    public static void moreGenePairs(String fileGenePairs,String fileGeneExpr,String filemiRExpr){
        ArrayList<String[]> listPairs=new ArrayList<>();
        String pairSplit="\t";
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


        ReadExpressionData geneExpr=new ReadExpressionData(geneNames);
        geneExpr.skipFirst=true;
        geneExpr.numberOfSamples =20531;
        geneExpr.separator="\t";
        geneExpr.genesAreRows=true;
        geneExpr.modifyName=false;
        geneExpr.readFile(fileGeneExpr);

        ArrayList<String> namesMiRNA= ReadFiles.extractNamesFromLine(filemiRExpr," ",0);
        ReadExpressionData miRExpr=new ReadExpressionData(namesMiRNA);
        miRExpr.skipFirst=true;
        miRExpr.separator=" ";
        miRExpr.genesAreRows=false;
        miRExpr.modifyName=false;
        miRExpr.readFile(filemiRExpr);

        int numberOfPermutations=1000;

           long timeStart=System.currentTimeMillis();
        for (int i = 0; i < listPairs.size(); i++) {
            String gene1=listPairs.get(i)[0];
            String gene2=listPairs.get(i)[1];
            String outputFileName="data/outputCMI"+gene1+"-"+gene2+".txt";
            computeCMIAndStore(geneExpr,miRExpr,gene1,gene2,outputFileName,numberOfPermutations);
        }

        long end = System.currentTimeMillis();
        System.out.println("time: " + (end - timeStart));


    }


    /**
     * Given two instances of ReadExpressionData, one with gene expression, the other with miRNA expression, computes
     * CMI values and stores them.
     * @param geneData ReadExpressionData instance with gene expression data stored
     * @param miRNAData ReadExpressionData instance with miRNA expression data stored
     * @param gene1 name of first gene
     * @param gene2 name of second gene
     * @param outputFileName name of output file
     * @param numberOfPermutations number of permutations to use for p-value
     */
    public static void computeCMIAndStore(ReadExpressionData geneData,ReadExpressionData miRNAData,String gene1,String gene2,String outputFileName,int numberOfPermutations){
        int gene1Index=geneData.getNameToData().get(gene1);
        int gene2Index=geneData.getNameToData().get(gene2);
        computeCMIAndStore(geneData,miRNAData,gene1Index,gene2Index,outputFileName,numberOfPermutations);
    }


    /**
     Given two instances of ReadExpressionData, one with gene expressioin, the other with miRNA expression, computes
     * CMI values and stores them.
     * @param geneData ReadExpressionData instance with gene expression data stored
     * @param miRNAData ReadExpressionData instance with miRNA expression data stored
     * @param gene1 index of first gene in geneData
     * @param gene2 index of second gene in geneData
     * @param outputFileName name of output file
     * @param numberOfPermutations number of permutations to use for p-value
     */
    public static void computeCMIAndStore(ReadExpressionData geneData,ReadExpressionData miRNAData,int gene1,int gene2,String outputFileName,int numberOfPermutations){
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
                ArrayList<double[]> data = new ArrayList<>();
                data.add(gene1Data);
                data.add(exprDatamiRna.get(i));
                data.add(gene2Data);
                CMIComplete cmiComplete = new CMIComplete(numberOfPermutations, data);
                //cmiComplete.cmiAndPValuePseudoUniform(8);
                cmiComplete.computeCMIandPValueBetterPartitioning();
                bw.write(miRNAData.getIntegersToNames().get(i)+"\t"+cmiComplete.cmi+"\t"+cmiComplete.pValue+"\n");
                bw.flush();
            }
            bw.write("\nI("+gene2Name+",miRNA|"+gene1Name+")\n");
            bw.write("miRNA id\tCMI\tp-value\n");
            for (int i = 0; i < exprDatamiRna.size(); i++) {
                ArrayList<double[]> data = new ArrayList<>();
                data.add(gene2Data);
                data.add(exprDatamiRna.get(i));
                data.add(gene1Data);
                CMIComplete cmiComplete = new CMIComplete(numberOfPermutations, data);
                //cmiComplete.cmiAndPValuePseudoUniform(8);
                cmiComplete.computeCMIandPValueBetterPartitioning();
                bw.write(miRNAData.getIntegersToNames().get(i)+"\t"+cmiComplete.cmi+"\t"+cmiComplete.pValue+"\n");
                bw.flush();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }




}
