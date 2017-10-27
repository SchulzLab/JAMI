package org.mpii.jami;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.mpii.jami.helpers.SettingsManager;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by fuksova on 10/7/15.
 */
public class Main {

    private static final Logger logger = LogManager.getLogger("JAMI");

    @Option(name="-output",usage="output file",metaVar="OUTPUT")
    File outputFile = new File("JAMI_CMI_results.txt");

    @Option(name="-set",usage="set if set notation should be used " +
            "as opposed to defining individual triplets to be tested",metaVar="INPUT")
    boolean setFormat;

    @Argument(index = 2,usage="file defining possible interactions between genes and miRNAs (set format use -set) or " +
            "triplets of gene-gene-miRNA",metaVar="INPUT")
    File genesMiRNA;

    @Argument(index = 1,usage="miRNA expression data",metaVar="INPUT")
    File filemiRNAExpr;

    @Argument(index = 0,usage="gene expression data",metaVar="INPUT")
    File fileGeneExpr;

    @Option(name="-method", hidden = true, usage = "for using an alternative CMI method. One of cupid, uniform, pseudouniform.")
    String method = "";

    @Option(name="-bins", hidden = true, usage = "number of bins when using uniform or pseudouniform CMI computation.")
    int numberOfBins = 0;

    @Option(name="-perm",usage="number of permutations for inferring empirical p-values. defaults to 1000.",metaVar="OPTIONS")
    int numberOfPermutations = 1000;

    @Option(name="-threads",usage="number of threads to use. -1 to use one less than the number of  available CPU cores",metaVar="OPTIONS")
    int numberOfThreads = -1;

    @Option(name="-noheader",usage="set this option if the input expression files have no headers",metaVar="OPTIONS")
    boolean noheader = false;

    @Option(name="-pcut",usage="optional p-value cutoff")
    double pValueCutoff = 1;

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

            if(fileGeneExpr == null || !fileGeneExpr.exists()){
                throw new IOException("ERROR: no gene expression file provided");
            }

            if(filemiRNAExpr == null || !filemiRNAExpr.exists()){
                throw new IOException("ERROR: no miRNA expression file provided");
            }

            if(genesMiRNA == null || !genesMiRNA.exists()){
                throw new IOException("ERROR: no gene expression file provided");
            }

            if(pValueCutoff <= 0 | pValueCutoff >= 1){
                throw new IllegalArgumentException("p-value cutoff -pcut has to be between 0 and 1.");
            }

            HashMap<String, Object> settings = new HashMap<>();
            settings.put("header", !noheader);
            settings.put("method", method);
            settings.put("numberOfBins", numberOfBins);
            settings.put("tripleFormat", !setFormat);
            settings.put("numberOfThreads", numberOfThreads);
            settings.put("numberOfPermutations", numberOfPermutations);
            settings.put("pValueCutoff", pValueCutoff);
            SettingsManager settingsManager = new SettingsManager(settings);

            CompleteRun completeRun = new CompleteRun(genesMiRNA,fileGeneExpr,filemiRNAExpr,outputFile, settingsManager);
            completeRun.runComputation();

        } catch(Exception e){
            System.err.println(e.getMessage());
            System.err.println("java JAMI [options...] gene_expression_file mir_expression_file gene_mir_interactions");
            parser.printUsage(System.err);
            System.err.println();

            return;
        }
    }
}


