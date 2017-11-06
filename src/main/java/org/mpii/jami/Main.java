package org.mpii.jami;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.spi.StringArrayOptionHandler;
import org.mpii.jami.helpers.SettingsManager;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * Created by fuksova on 10/7/15.
 */
public class Main {

    private static final Logger logger = LogManager.getLogger("JAMI");

    @Option(name="-output",usage="output file")
    File outputFile = new File("JAMI_CMI_results.txt");

    @Option(name="-set",usage="set if set notation should be used " +
            "as opposed to defining individual triplets to be tested")
    boolean setFormat;

    @Argument(index = 2,usage="file defining possible interactions between genes and miRNAs (set format use -set) or " +
            "triplets of gene-gene-miRNA")
    File genesMiRNA;

    @Argument(index = 1,usage="miRNA expression data")
    File filemiRNAExpr;

    @Argument(index = 0,usage="gene expression data")
    File fileGeneExpr;

    @Option(name="-method", hidden = true, usage = "for using an alternative CMI method. One of cupid, uniform, pseudouniform.")
    String method = "";

    @Option(name="-bins", hidden = true, usage = "number of bins when using uniform or pseudouniform CMI computation.")
    int numberOfBins = 0;

    @Option(name="-perm",usage="number of permutations for inferring empirical p-values. defaults to 1000.")
    int numberOfPermutations = 1000;

    @Option(name="-threads",usage="number of threads to use. -1 to use one less than the number of  available CPU cores")
    int numberOfThreads = -1;

    @Option(name="-noheader",usage="set this option if the input expression files have no headers")
    boolean noheader = false;

    @Option(name="-pcut",usage="optional Benjamini Hochberg adjusted p-value cutoff")
    double pValueCutoff = 1.0;

    @Option(name="-genes", handler = StringArrayOptionHandler.class, usage="filter for miRNA triplets with this gene as regulator")
    List<String> genes;

    @Option(name="-restricted",usage="set this option to restrict analysis to interactions between the selected gene")
    boolean restricted = false;


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

            if(parser.getArguments().size() != 3){
                throw new IOException("ERROR: number of arguments does not match");
            }

            if(fileGeneExpr == null || !fileGeneExpr.exists()){
                throw new IOException("ERROR: no gene expression file provided");
            }

            if(filemiRNAExpr == null || !filemiRNAExpr.exists()){
                throw new IOException("ERROR: no miRNA expression file provided");
            }

            if(genesMiRNA == null || !genesMiRNA.exists()){
                throw new IOException("ERROR: no gene expression file provided");
            }

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
            settings.put("selectedGenes", new HashSet<>(genes));
            settings.put("restricted", restricted);
            SettingsManager settingsManager = new SettingsManager(settings);

            CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,outputFile, settingsManager);
            completeRun.runComputation();

        } catch(Exception e){
            System.err.println("ERROR:" + e.getMessage());
            System.err.println("ERROR DETAILS:");
            e.printStackTrace(System.err);
            System.err.println("JAMI USAGE:");
            System.err.println("java JAMI [options...] gene_expression_file mir_expression_file gene_mir_interactions");
            parser.printUsage(System.err);
            System.err.println();

            return;
        }
    }
}


