package org.mpii.jami;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.spi.StringArrayOptionHandler;
import org.mpii.jami.helpers.HelpOptionHandler;
import org.mpii.jami.helpers.SettingsManager;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;

/**
 * Created by fuksova on 10/7/15.
 */
public class Main {

    private static final Logger logger = LogManager.getLogger("JAMI");

    @Option(name="-output",usage="output file")
    File outputFile = new File("JAMI_CMI_results.txt");

    @Option(name="-set",usage="set if set notation should be used " +
            "as opposed to defining individual triplets to be tested",
            handler = HelpOptionHandler.class)
    boolean setFormat;

    @Argument(index = 2,usage="file defining possible interactions between genes and miRNAs (set format use -set) or " +
            "triplets of gene-gene-miRNA")
    File genesMiRNA;

    @Argument(index = 1,usage="miRNA expression data")
    File filemiRNAExpr;

    @Argument(index = 0, usage="gene expression data")
    File fileGeneExpr;

    @Option(name="-method", hidden = true, usage = "for using an alternative CMI method. One of cupid, uniform, pseudouniform.")
    String method = SettingsManager.defaultMethod;

    @Option(name="-nozeros", usage = "set this flag to ignore duplicated zero expression values. Otherwise JAMI" +
            "will group duplicates of the lowest values such as 0 or the smallest negative value in log scaled data.")
    boolean nozeros = false;

    @Option(name="-bins", hidden = true, usage = "number of bins when using uniform or pseudouniform CMI computation.")
    int numberOfBins = 0;

    @Option(name="-perm",usage="number of permutations for inferring empirical p-values.")
    int numberOfPermutations = SettingsManager.defaultPermutations;

    @Option(name="-threads",usage="number of threads to use. -1 to use one less than the number of available CPU cores")
    int numberOfThreads = SettingsManager.defaultThreads;

    @Option(name="-batch",usage="number of triplets in each batch. affects overhead of multi-threaded computation")
    int batchSize = SettingsManager.defaultBatchSize;

    @Option(name="-noheader",usage="set this option if the input expression files have no headers",
            handler = HelpOptionHandler.class)
    boolean noheader = false;

    @Option(name="-pcut",usage="optional Benjamini Hochberg adjusted p-value cutoff")
    double pValueCutoff = SettingsManager.defaultPValueCutoff;

    @Option(name="-pchi", usage = "significance level for the chi-squared test in adaptive partitioning")
    double pChiSquare = SettingsManager.defaultPChiSquare;

    @Option(name="-genes", handler = StringArrayOptionHandler.class, usage="filter for miRNA triplets with this gene or these genes as regulator")
    List<String> genes;

    @Option(name="-restricted",usage="set this option to restrict analysis to interactions between the selected genes",
            handler = HelpOptionHandler.class)
    boolean restricted = false;

    @Option(name="-v", usage="show JAMI version",
            handler = HelpOptionHandler.class)
    boolean showVersion;

    @Option(name="-h", usage="show this usage information",
            handler = HelpOptionHandler.class)
    boolean showHelp;

    @Option(name="-verbose", usage="show verbose error messages",
            handler = HelpOptionHandler.class)
    boolean verboseLogging = false;

    /**
     * Uncomment one of the options
     * @param args
     */
    public static void main(String []args) throws IOException {

        logger.info("Java Conditional Mutual Information (JAMI) started.");

        new Main().runJAMI(args);
        logger.info("Exit");
    }

    public void runJAMI(String[] args) throws IOException{
        CmdLineParser parser = new CmdLineParser(this);

        try {
            parser.parseArgument(args);

            if(showHelp){
                System.out.println("You may find a detailed documentation at http://jami.readthedocs.io/en/latest/index.html");
                System.out.println("JAMI USAGE:");
                System.out.println("java JAMI [options...] gene_expression_file mir_expression_file gene_mir_interactions");
                parser.printUsage(System.out);
                System.exit(0);
            }

            if (showVersion) {
                Properties prop = new Properties();
                ClassLoader classLoader = getClass().getClassLoader();
                prop.load(classLoader.getResourceAsStream("config.properties"));
                System.out.println("JAMI " + prop.getProperty("version") + " (BUILD " + prop.getProperty("build") + ")");
                System.exit(0);
            }

            boolean inputOk = true;
            if(fileGeneExpr == null || !fileGeneExpr.exists()){
                logger.error("no gene expression file provided");
                inputOk = false;
            }

            if(filemiRNAExpr == null || !filemiRNAExpr.exists()){
                logger.error("no miRNA expression file provided");
                inputOk = false;
            }

            if(genesMiRNA == null || !genesMiRNA.exists()){
                logger.error("no gene-miRNA interactions file provided");
                inputOk = false;
            }

            if(!inputOk) throw new IllegalArgumentException("Arguments missing or incorrect.");

            if(pValueCutoff <= 0.0 && pValueCutoff >= 1.0){
                throw new IllegalArgumentException("p-value cutoff -pcut has to be between 0 and 1.");
            }

            if(!setFormat && restricted){
                throw new IllegalArgumentException("restricted option is only valid for the set input format.");
            }

            HashMap<String, Object> settings = new HashMap<>();
            settings.put("header", !noheader);
            settings.put("method", method);
            settings.put("numberOfBins", numberOfBins);
            settings.put("tripleFormat", !setFormat);
            settings.put("numberOfThreads", numberOfThreads);
            settings.put("numberOfPermutations", numberOfPermutations);
            settings.put("pValueCutoff", pValueCutoff);
            settings.put("batchSize", batchSize);
            settings.put("considerZeros", !nozeros);
            settings.put("pChiSquare", pChiSquare);
            if(genes != null) settings.put("selectedGenes", new HashSet<>(genes));
            else settings.put("selectedGenes", null);
            settings.put("restricted", restricted);
            SettingsManager settingsManager = new SettingsManager(settings);

            CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,outputFile, settingsManager);
            completeRun.runComputation();

        } catch(Exception e){
            logger.error(e.getMessage());
            if(verboseLogging) {
                logger.error("DETAILS:", e);
            }
            System.out.println("JAMI USAGE:");
            System.out.println("java -jar JAMI [options...] gene_expression_file mir_expression_file gene_mir_interactions [options...]");
            parser.printUsage(System.out);
            System.out.println();

            return;
        }
    }
}


