package org.mpii.jami;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ForkJoinPool;

/**
 * Created by fuksova on 3/14/17.
 */
public class CompleteRun {
    private static final Logger logger = LogManager.getLogger("JAMI");
    private Character separator = '\t';
    public boolean completed;
    public int tripletsWrittenToDisk;
    private File outputFile;
    private File fileGenesMiRNA;
    private File fileGeneExpr;
    private File filemiRExpr;
    private int numberOfPermutations;
    private ExpressionData geneExpr;
    private ExpressionData miRExpr;
    private BufferedWriter bw;
    private boolean tripleFormat;
    private String method = "";
    private int numberOfBins = 0;
    private int numberOfThreads = -1;

    public CompleteRun(File fileGenesMiRNA,File fileGeneExpr,File filemiRExpr,File outputFile,
                       int numberOfPermutations,boolean tripleFormat,
                       String method, int numberOfBins, int numberOfThreads){
        this.tripletsWrittenToDisk = 0;
        this.completed = false;
        this.method = method;
        this.numberOfBins = numberOfBins;
        this.outputFile=outputFile;
        this.fileGenesMiRNA=fileGenesMiRNA;
        this.fileGeneExpr=fileGeneExpr;
        this.filemiRExpr=filemiRExpr;
        this.tripleFormat=tripleFormat;
        this.numberOfThreads=numberOfThreads;
        this.numberOfPermutations=numberOfPermutations;
    }

    public CompleteRun(File fileGenesMiRNA,File fileGeneExpr,File filemiRExpr,File outputFile,
                       int numberOfPermutations,boolean tripleFormat){
        this.tripletsWrittenToDisk = 0;
        this.completed = false;
        this.outputFile=outputFile;
        this.fileGenesMiRNA=fileGenesMiRNA;
        this.fileGeneExpr=fileGeneExpr;
        this.filemiRExpr=filemiRExpr;
        this.tripleFormat=tripleFormat;
        this.numberOfPermutations=numberOfPermutations;
    }


    /**
     * Reads the files and performs all CMI computations.
     */
    public void runComputation(){

        if(this.completed){
            logger.error("Computation on this object was previously finished.");
        }

        InteractionData interactions = new InteractionData();

        if(tripleFormat){
            logger.debug("Reading CMI candidate file in triplet format.");

            interactions.readFileWithTriples(this.fileGenesMiRNA);
        }
        else{
            logger.debug("Reading CMI candidate file in set format.");

            interactions.readFileInSetFormat(this.fileGenesMiRNA);
        }

        //read only gene and miRNA expression data we actually need
        ArrayList<String> geneNames=new ArrayList<>(interactions.getGenes());
        ArrayList<String> miRNANames=new ArrayList<>(interactions.getMiRNAs());

        miRExpr=new ExpressionData(miRNANames);
        geneExpr=new ExpressionData(geneNames);

        geneExpr.readFile(fileGeneExpr);
        miRExpr.readFile(filemiRExpr);

        CMIComplete.initRandomized(geneExpr.getNumberOfSamples(), numberOfPermutations);

        //parallelization
        ForkJoinPool fjpool;

        if(numberOfThreads == -1)
             fjpool = ForkJoinPool.commonPool();
        else fjpool = new ForkJoinPool(numberOfThreads);

        long timeStart=System.currentTimeMillis();

        FileWriter fw;
        try {
            fw = new FileWriter(outputFile);
            bw = new BufferedWriter(fw);
            bw.write("Interacting");
            bw.write(separator);
            bw.write("Conditional");
            bw.write(separator);
            bw.write("Mediator");
            bw.write(separator);
            bw.write("CMI");
            bw.write(separator);
            bw.write("p-value\n");

            for (String[] oneTriple : interactions.getTriplets()) {
                Integer gene1Index = geneExpr.getNameToData().get(oneTriple[0]);
                Integer gene2Index = geneExpr.getNameToData().get(oneTriple[1]);
                Integer miRNAIndex = miRExpr.getNameToData().get(oneTriple[2]);

                if (gene1Index != null && gene2Index != null && miRNAIndex != null) {
                    computeCMIAndStore(gene1Index, gene2Index, miRNAIndex, fjpool);
                    this.tripletsWrittenToDisk++;
                }

                if(gene1Index == null){
                    logger.debug("Gene " + oneTriple[0] + " not found in gene expression data");
                }
                if(gene2Index == null){
                    logger.debug("Gene " + oneTriple[1] + " not found in gene expression data");
                }
                if(miRNAIndex == null){
                    logger.debug("miRNA " + oneTriple[2] + " not found in miRNA expression data");
                }

            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        this.completed = true;


        long end = System.currentTimeMillis();
        logger.info("Computation finished in " + (end - timeStart)/1000 + "s");


    }


    /**
     * Computes CMI of genes and miRNA. Both (gene1,miRNA|gene2) and (gene2,miRNA|gene1) are computed and written into
     * the output file.
     * @param gene1 index of gene1 in geneExpr
     * @param gene2 index of gene2 in geneExpr
     * @param miRNA index of miRNA in miExpr
     */
    public void computeCMIAndStore(int gene1, int gene2, int miRNA,
                                   ForkJoinPool fjpool) {
        double[] gene1Data = geneExpr.getExpressionData().get(gene1);
        double[] gene2Data = geneExpr.getExpressionData().get(gene2);
        double[] miRNAData=miRExpr.getExpressionData().get(miRNA);
        String gene1Name = geneExpr.getIntegersToNames().get(gene1);
        String gene2Name = geneExpr.getIntegersToNames().get(gene2);
        String miRNAName = miRExpr.getIntegersToNames().get(miRNA);


        try {

            bw.write(gene1Name + separator + gene2Name + separator + miRNAName + separator);
            double[] result = computeTriple(gene2Data, gene1Data, miRNAData, fjpool);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();

            bw.write(gene2Name + separator + gene1Name + separator + miRNAName + separator);
            result = computeTriple(gene1Data, gene2Data, miRNAData, fjpool);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();


        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    /**
     * Run computation of CMI (geneInterData,miRNAData|geneCondData).
     * @param geneCondData expression data
     * @param geneInterData expression data
     * @param miRNAData expression data
     * @return First element is CMI, second is p-value;
     */
    private double[] computeTriple(double[] geneCondData, double[] geneInterData, double[] miRNAData,
                                   ForkJoinPool fjpool) {  //gene1Data is randomized
        ArrayList<double[]> data = new ArrayList<>();
        double[] result=new double[2];
        data.add(geneInterData);
        data.add(miRNAData);
        data.add(geneCondData);  //this is randomized
        CMIComplete cmiComplete;

        cmiComplete = new CMIComplete(numberOfPermutations, data, fjpool);

        switch(method){
            case "cupid": cmiComplete.computeAsCUPID();
                break;
            case "uniform": cmiComplete.computeUniformGrid(numberOfBins);
                break;
            case "pseudouniform": cmiComplete.computePseudoUniformGrid(numberOfBins);
                break;
            default: cmiComplete.computeIterativePartitioning();
        }
        result[0]=cmiComplete.cmi;
        result[1]=cmiComplete.pValue;
        return result;

    }




}
