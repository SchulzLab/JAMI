package org.mpii.jami;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mpii.jami.model.Triplet;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ForkJoinPool;

/**
 * Created by fuksova on 3/14/17.
 */
public class CompleteRun{
    private static final Logger logger = LogManager.getLogger("JAMI");
    private String separator = "\t";
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
    private boolean header;
    private HashMap<Triplet, Double> cmis = new HashMap<>();
    private HashMap<Triplet, Double> pvalues = new HashMap<>();

    public CompleteRun(File fileGenesMiRNA,File fileGeneExpr,File filemiRExpr,File outputFile,
                       int numberOfPermutations,boolean tripleFormat,
                       String method, int numberOfBins, int numberOfThreads,
                       boolean header){
        this.tripletsWrittenToDisk = 0;
        this.completed = false;
        this.header = header;
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
                       int numberOfPermutations,boolean tripleFormat,
                       boolean header){
        this.tripletsWrittenToDisk = 0;
        this.completed = false;
        this.outputFile=outputFile;
        this.fileGenesMiRNA=fileGenesMiRNA;
        this.fileGeneExpr=fileGeneExpr;
        this.filemiRExpr=filemiRExpr;
        this.tripleFormat=tripleFormat;
        this.numberOfPermutations=numberOfPermutations;
        this.header = header;
    }

    public HashMap<Triplet, Double> getCmis() {
        return cmis;
    }

    public HashMap<Triplet, Double> getPvalues() {
        return pvalues;
    }

    /**
     * Reads the files and performs all CMI computations.
     */
    public void runComputation() throws IOException {

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

        geneExpr.readFile(fileGeneExpr, header);
        miRExpr.readFile(filemiRExpr, header);

        if(geneExpr.getNumberOfSamples() != miRExpr.getNumberOfSamples())
            throw new IOException("Gene and miRNA expression files have differing sample numbers.");

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
            bw.write("Modulator");
            bw.write(separator);
            bw.write("Target");
            bw.write(separator);
            bw.write("miRNA");
            bw.write(separator);
            bw.write("CMI");
            bw.write(separator);
            bw.write("p-value\n");

            for (Triplet t : interactions.getTriplets()) {
                computeCMIAndStore(t, fjpool);
                this.tripletsWrittenToDisk+=2;
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
     * @param t The triplet for which we want to compute the CMI
     * @param fjpool ForkJoinPool to use for parallel execution
     */
    public void computeCMIAndStore(Triplet t, ForkJoinPool fjpool) {

        String gene1Name = t.getGeneOne();
        String gene2Name = t.getGeneTwo();
        String miRNAName = t.getMiRNA();

        Integer gene1Index = geneExpr.getNameToData().get(gene1Name);
        Integer gene2Index = geneExpr.getNameToData().get(gene2Name);
        Integer miRNAIndex = miRExpr.getNameToData().get(miRNAName);

        if(gene1Index == null){
            logger.warn("Gene " + gene1Name + " not found in gene expression data");
            return;
        }
        if(gene2Index == null){
            logger.warn("Gene " + gene2Name + " not found in gene expression data");
            return;
        }
        if(miRNAIndex == null){
            logger.warn("miRNA " + miRNAName + " not found in miRNA expression data");
            return;
        }

        double[] gene1Data = geneExpr.getExpressionData().get(gene1Index);
        double[] gene2Data = geneExpr.getExpressionData().get(gene2Index);
        double[] miRNAData=miRExpr.getExpressionData().get(miRNAIndex);

        try {
            bw.write(gene1Name + separator + gene2Name + separator + miRNAName + separator);
            double[] result = computeTriple(gene1Data, gene2Data, miRNAData, fjpool);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();

            this.cmis.put(t, result[0]);
            this.pvalues.put(t, result[1]);

            bw.write(gene2Name + separator + gene1Name + separator + miRNAName + separator);
            result = computeTriple(gene2Data, gene1Data, miRNAData, fjpool);
            bw.write(result[0] + separator + result[1] + "\n");
            bw.flush();

            this.cmis.put(new Triplet(gene2Name, gene1Name, miRNAName), result[0]);
            this.pvalues.put(new Triplet(gene2Name, gene1Name, miRNAName), result[1]);

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
