package org.mpii.jami;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mpii.jami.cmi.Cube;
import org.mpii.jami.cmi.IterativePartitioning;
import org.mpii.jami.helpers.*;
import org.mpii.jami.input.ExpressionData;
import org.mpii.jami.input.InteractionData;
import org.mpii.jami.model.Triplet;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by fuksova on 3/14/17. modified by mlist
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
    private double pValueCutoff = 1;
    private HashSet<String> selectedGenes;
    private List<Triplet> triplets;
    private List<Triplet> omittedTriplets;
    private Cube initialCube;
    final int batchSize = 100;
    final double progressMinDiff = 0.01;
    private boolean restricted = false;
    private NonRepetitiveLogger nonRepetitiveLogger = new NonRepetitiveLogger();

    public CompleteRun(File fileGenesMiRNA, File fileGeneExpr, File filemiRExpr, File outputFile,
                       SettingsManager settingsManager){
        this.tripletsWrittenToDisk = 0;
        this.completed = false;
        this.header = (boolean) settingsManager.get("header");
        this.method = (String) settingsManager.get("method");
        this.numberOfBins = (int) settingsManager.get("numberOfBins");
        this.outputFile = outputFile;
        this.fileGenesMiRNA = fileGenesMiRNA;
        this.fileGeneExpr = fileGeneExpr;
        this.filemiRExpr = filemiRExpr;
        this.tripleFormat = (boolean) settingsManager.get("tripleFormat");
        this.numberOfThreads = (int) settingsManager.get("numberOfThreads");
        this.numberOfPermutations = (int) settingsManager.get("numberOfPermutations");
        this.pValueCutoff = (double) settingsManager.get("pValueCutoff");
        this.selectedGenes = (HashSet<String>) settingsManager.get("selectedGenes");
        this.restricted = (boolean) settingsManager.get("restricted");
    }

    public List<Triplet> getTriplets() {
        return this.triplets;
    }

    public List<Triplet> getOmittedTriplets() {
        return omittedTriplets;
    }

    public CompleteRun(File fileGenesMiRNA, File fileGeneExpr, File filemiRExpr, File outputFile){
        this(fileGenesMiRNA, fileGeneExpr, filemiRExpr, outputFile, new SettingsManager());
    }

    public HashSet<String> getSelectedGenes() {
        return this.selectedGenes;
    }

    public void filterForGene(String selectedGene) {
        this.selectedGenes = new HashSet<>();
        this.selectedGenes.add(selectedGene);
    }

    /**
     * Reads the files and performs all CMI computations.
     */
    public void runComputation() throws IOException, ExecutionException, InterruptedException {

        if(this.completed){
            logger.error("Computation on this object was previously finished.");
        }

        InteractionData interactions = new InteractionData();

        if(tripleFormat){
            logger.info("Reading CMI candidate file in triplet format.");

            interactions.readFileWithTriples(this.fileGenesMiRNA, selectedGenes);
        }
        else{
            logger.info("Reading CMI candidate file in set format.");

            interactions.readFileInSetFormat(this.fileGenesMiRNA, selectedGenes, restricted);
        }

        logger.info("" + interactions.getTriplets().size() + " interactions (triplets) selected for processing.");

        if(interactions.getTriplets().size() == 0){
            logger.info("Exit");
            System.exit(0);
        }

        //read only gene and miRNA expression data we actually need
        ArrayList<String> geneNames=new ArrayList<>(interactions.getGenes());
        ArrayList<String> miRNANames=new ArrayList<>(interactions.getMiRNAs());

        miRExpr=new ExpressionData(miRNANames);
        geneExpr=new ExpressionData(geneNames);

        logger.info("Reading gene expression data from " + fileGeneExpr);
        geneExpr.readFile(fileGeneExpr, header);

        logger.info("Reading miRNA expression data from " + filemiRExpr);
        miRExpr.readFile(filemiRExpr, header);

        if(geneExpr.getNumberOfSamples() != miRExpr.getNumberOfSamples())
            throw new IOException("Gene and miRNA expression files have differing sample numbers.");

        initialCube = IterativePartitioning.getInitialCube(geneExpr.getNumberOfSamples());

        Map<Boolean, List<Triplet>> indexedTriplets = interactions.getTriplets().stream().map(t -> {
            String gene1Name = t.getGeneOne();
            String gene2Name = t.getGeneTwo();
            String miRNAName = t.getMiRNA();

            Integer gene1Index = geneExpr.getNameToData().get(gene1Name);
            Integer gene2Index = geneExpr.getNameToData().get(gene2Name);
            Integer miRNAIndex = miRExpr.getNameToData().get(miRNAName);

            if(gene1Index == null){
                nonRepetitiveLogger.warn("Gene " + gene1Name + " not found in gene expression data");
                t.markInvalid();
                return t;
            }
            if(gene2Index == null){
                nonRepetitiveLogger.warn("Gene " + gene2Name + " not found in gene expression data");
                t.markInvalid();
                return t;
            }
            if(miRNAIndex == null){
                nonRepetitiveLogger.warn("miRNA " + miRNAName + " not found in miRNA expression data");
                t.markInvalid();
                return t;
            }

            t.setGeneOneIndex(gene1Index);
            t.setGeneTwoIndex(gene2Index);
            t.setMiRNAIndex(miRNAIndex);
            return(t);

        }).collect(Collectors.partitioningBy(Triplet.check()));

        this.triplets = indexedTriplets.get(true);
        this.omittedTriplets = indexedTriplets.get(false);

        logger.info("" + interactions.getTriplets().size() + " interactions (triplets) mapped to expression data.");

        CMIComplete.initRandomized(geneExpr.getNumberOfSamples(), numberOfPermutations);

        //parallelization
        ForkJoinPool fjpool;

        if(numberOfThreads == -1)
             fjpool = ForkJoinPool.commonPool();
        else fjpool = new ForkJoinPool(numberOfThreads);

        long timeStart=System.currentTimeMillis();

        logger.info("Computing CMI and p-values...");

        //the following scheduling is only for the progress bar
        //with this complicated construct we make sure that only the main thread is tasked with keeping track
        int numberOfElementsPerBatch = batchSize;
        int numberOfBatches = (this.triplets.size()+numberOfElementsPerBatch-1)/numberOfElementsPerBatch;
        AtomicInteger currentBatch = new AtomicInteger(0);
        ScheduledExecutorService es = Executors.newScheduledThreadPool(1);

        ScheduledFuture<?> f = es.scheduleWithFixedDelay(new Runnable() {
            double progress = 0;
            public void run() {
                double newProgress = (double) currentBatch.get() / numberOfBatches;

                if ((newProgress - progress) > progressMinDiff) {
                    progress = newProgress;
                    Progressbar.updateProgress(progress);
                }
            }
        }, 100, 100, TimeUnit.MILLISECONDS);

        List<List<Double>> immutableGeneExpression = Collections.unmodifiableList(geneExpr.getExpressionData());
        List<List<Double>> immutableMiRNAExpression = Collections.unmodifiableList(miRExpr.getExpressionData());

        List<List<Triplet>> batches = fjpool.submit(() ->
            IntStream.range(0, numberOfBatches)
                    .mapToObj(i -> this.triplets.subList(i*numberOfElementsPerBatch,
                            Math.min(this.triplets.size(), (i+1)*numberOfElementsPerBatch)))
                    .parallel()
                    .map(batch -> batch.stream()
                            .map(triplet -> computeCMI(triplet, immutableGeneExpression, immutableMiRNAExpression))
                            .collect(Collectors.toList()))
                    .peek(batch -> currentBatch.getAndIncrement())
                    .collect(Collectors.toList())
                    ).get();

        this.triplets = BenjaminiHochberg.adjustPValues(batches.stream()
                .flatMap(batch -> batch.stream())
                .collect(Collectors.toList()));

        f.cancel(true);
        es.shutdown();

        this.completed = true;

        long end = System.currentTimeMillis();
        logger.info("Computation finished in " + (end - timeStart)/1000 + "s");

        logger.info("Saving results to " + outputFile);
        saveResults(this.triplets);
    }

    private void saveResults(List<Triplet> triplets) {
        logger.debug("Writing triplets to output file");
        FileWriter fw;
        try {
            fw = new FileWriter(outputFile);
            bw = new BufferedWriter(fw);
            bw.write("Source");
            bw.write(separator);
            bw.write("Target");
            bw.write(separator);
            bw.write("miRNA");
            bw.write(separator);
            bw.write("CMI");
            bw.write(separator);
            bw.write("p-value");
            bw.write(separator);
            bw.write("p-adjusted\n");

            for(Triplet t : triplets){
                if(t.getpAdjust() <= this.pValueCutoff) {
                    bw.write(t.getGeneOne() + separator + t.getGeneTwo() + separator + t.getMiRNA() + separator);
                    bw.write(t.getCmi() + separator + t.getpValue() + separator + t.getpAdjust() + "\n");
                    this.tripletsWrittenToDisk++;
                }
            }

            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        logger.debug("Writing aggregated gene pairs to output file");

        try {
            String parent = "";
            if(outputFile.getParent() != null) parent = outputFile.getParent();
            fw = new FileWriter(new File(parent + "aggregated_" + outputFile.getName()));
            bw = new BufferedWriter(fw);
            bw.write("Source");
            bw.write(separator);
            bw.write("Target");
            bw.write(separator);
            bw.write("Number of miRNAs");
            bw.write(separator);
            bw.write("miRNAs");
            bw.write(separator);
            bw.write("CMI min");
            bw.write(separator);
            bw.write("CMI max");
            bw.write(separator);
            bw.write("p-value min");
            bw.write(separator);
            bw.write("p-value min");
            bw.write(separator);
            bw.write("p-adjusted min");
            bw.write(separator);
            bw.write("p-adjusted max");
            bw.write(separator);
            bw.write("p-adjusted Fisher\n");

            triplets.stream()
                    .collect(Collectors.groupingBy(Triplet::getGeneOne,
                            Collectors.groupingBy(Triplet::getGeneTwo,
                                    Collectors.toList())))
                    .forEach((geneOne, stringListMap) -> {
                        stringListMap
                                .forEach((geneTwo, listOfTriplets) -> {
                                    try {
                                        double pFisher = FisherMethod.combinedPValue(
                                                listOfTriplets.stream()
                                                        .map(t -> t.getpAdjust())
                                                        .collect(Collectors.toList()));
                                        if(pFisher > pValueCutoff) {
                                            bw.write(geneOne);
                                            bw.write(separator);
                                            bw.write(geneTwo);
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.size()));
                                            bw.write(separator);
                                            bw.write(listOfTriplets.stream()
                                                    .map(t -> t.getMiRNA())
                                                    .collect(Collectors.joining("/"))
                                            );
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getCmi())
                                                    .min().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getCmi())
                                                    .max().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getpValue())
                                                    .min().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getpValue())
                                                    .max().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getpAdjust())
                                                    .min().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(listOfTriplets.stream()
                                                    .mapToDouble(t -> t.getpAdjust())
                                                    .max().getAsDouble()));
                                            bw.write(separator);
                                            bw.write(String.valueOf(pFisher) + "\n");
                                        }
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                });
                    });
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * Computes CMI of genes and miRNA. Both (gene1,miRNA|gene2) and (gene2,miRNA|gene1) are computed and written into
     * the output file.
     * @param t The triplet for which we want to compute the CMI
     */
    public Triplet computeCMI(Triplet t, List<List<Double>> geneExpr, List<List<Double>> miRExpr) {

        List<Double> gene1Data = geneExpr.get(t.getGeneOneIndex());
        List<Double> gene2Data = geneExpr.get(t.getGeneTwoIndex());
        List<Double> miRNAData = miRExpr.get(t.getMiRNAIndex());

        double[] result = computeTriple(gene1Data, gene2Data, miRNAData);

        t.setCmi(result[0]);
        t.setpValue(result[1]);

        return(t);
    }

    /**
     * Run computation of CMI (geneInterData,miRNAData|geneCondData).
     * @param geneCondData expression data
     * @param geneInterData expression data
     * @param miRNAData expression data
     * @return First element is CMI, second is p-value;
     */
    private double[] computeTriple(List<Double> geneCondData, List<Double> geneInterData, List<Double> miRNAData) {  //gene1Data is randomized
        ArrayList<List<Double>> data = new ArrayList<>();
        double[] result=new double[2];
        data.add(geneInterData);
        data.add(miRNAData);
        data.add(geneCondData);  //this is randomized
        CMIComplete cmiComplete;

        cmiComplete = new CMIComplete(numberOfPermutations, data);

        String method = "default";
        switch(method){
            case "cupid": cmiComplete.computeAsCUPID();
                break;
            case "uniform": cmiComplete.computeUniformGrid(this.numberOfBins);
                break;
            case "pseudouniform": cmiComplete.computePseudoUniformGrid(this.numberOfBins);
                break;
            default: cmiComplete.computeIterativePartitioning(initialCube);
        }
        result[0]=cmiComplete.cmi;
        result[1]=cmiComplete.pValue;
        return result;

    }




}
