package org.mpii.jami.input;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.util.FileUtils;
import org.mpii.jami.model.Triplet;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import static java.util.Arrays.stream;

/**
 * Created by mlist on 10/24/17.
 */
public class InteractionData {

    private static final Logger logger = LogManager.getLogger("JAMI");

    private HashSet<String> genes = new HashSet<>();
    private HashSet<String> miRNAs = new HashSet<>();

    public HashSet<String> getGenes() {
        return genes;
    }

    public HashSet<String> getMiRNAs() {
        return miRNAs;
    }

    public HashSet<Triplet> getTriplets() { return triplets; }

    private HashSet<Triplet> triplets = new HashSet<>();

    public InteractionData(){

    }

    /**
     * Reads file with triple of two genes and one miRNA names whose interactions will be tested.
     * Every row contains one triple to test separated by "separator".
     * @param file Name of file
     * @param separator Column separator
     * @param skipFirst Skip first row
     * @return ArrayList with all triple from the file
     */
    private void readTriples(File file, String separator, boolean skipFirst, HashSet<String> selectedGenes){
        Pattern pattern = Pattern.compile(separator);
        HashMap<String, Boolean> selectedGenesFound = new HashMap<>();
        for(String gene : selectedGenes){
            selectedGenesFound.put(gene, false);
        }

        try {
            BufferedReader br = readFile(file);
            String line;
            if(skipFirst){
                br.readLine();
            }
            line=br.readLine();
            while (line!=null) {
                String[] entries = pattern.split(line);
                if((selectedGenes != null && selectedGenes.contains(entries[0])) || selectedGenes == null) {
                    Triplet entry = new Triplet(entries);
                    if (!this.triplets.contains(entry)) this.triplets.add(entry);
                    if(selectedGenes != null) selectedGenesFound.put(entries[0], true);
                    this.genes.add(entries[0]);
                    this.genes.add(entries[1]);
                    this.miRNAs.add(entries[2]);
                }
                if(selectedGenes == null) {
                    Triplet reverseEntry = new Triplet(entries[1], entries[0], entries[2]);
                    if (!this.triplets.contains(reverseEntry)) this.triplets.add(reverseEntry);
                }

                line=br.readLine();
            }
            br.close();
            selectedGenesFound
                    .keySet()
                    .stream()
                    .filter(gene -> !selectedGenesFound.get(gene))
                    .forEach(gene -> logger.warn("No miRNA interactions were found for " + gene));


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static BufferedReader readFile(File file) throws IOException {
        if(FileUtils.getFileExtension(file).equals("gz")) {
            FileInputStream fis = new FileInputStream(file);
            GZIPInputStream gis = new GZIPInputStream(fis);
            InputStreamReader isr = new InputStreamReader(gis);
            return(new BufferedReader(isr));
        } else{
            FileReader fr = new FileReader(file);
            return(new BufferedReader(fr));
        }
    }

    /**
     * Read files with genes-miRNA interactions. Every line contains one gene name and a set of miRNA names that
     * are assumed to interact with this gene. The format of one line is: geneID, "separatorMain" and list of
     * miRNA names separated by "separatorMiRNA".
     * @param file Name of input file.
     * @param separatorMain Separates the column with genes from the set of miRNA interacting with it.
     * @param separatorMiRNA Separates miRNA names in on set.
     * @param skipFirst Skip first line in the file.
     * @return HashMap geneID->array of miRNA IDs.
     */
    private static HashMap<String,String[]> geneToMiRNA(File file,String separatorMain,String separatorMiRNA,boolean skipFirst){
        Pattern pattern = Pattern.compile(separatorMain);
        Pattern patternMiRNA=Pattern.compile(separatorMiRNA);

        HashMap<String,String[]> genesWithMiRNA=new HashMap<>();
        try {

            BufferedReader br = readFile(file);

            String line=br.readLine();
            while (line!=null) {
                if(skipFirst){
                    line = br.readLine();
                    skipFirst = false;
                    continue;
                }
                String[] entries = pattern.split(line);

                if(entries.length == 1){
                    logger.debug("gene " + entries[0] + "has no miRNAs");
                }
                else if(entries.length!=2){  //maybe add possibility that there is no miRNA
                    throw new IOException("At least one line in file "+file +
                            " has more zero or more than two elements");
                }
                String geneName=entries[0];
                String[] allMiRNA=patternMiRNA.split(entries[1]);
                genesWithMiRNA.put(geneName,allMiRNA);
                line=br.readLine();
            }
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return genesWithMiRNA;
    }

    /**
     * Reads file with gene-miRNA interactions in the format where every gene has assigned a set of miRNA.
     * Then it creates triples of two genes and one miRNA to be tested based on common intersection of the miRNA sets
     * of corresponding to genes.
     * miRExpr and geneExpr are initialized with list of gene and miRNA names found in the file.
     */
    public void readFileInSetFormat(File fileGenesMiRNA, HashSet<String> selectedGenes, boolean restricted){
        HashMap<String, String[]> genesToMiRNA = geneToMiRNA(fileGenesMiRNA, "\t", ",", true);

        this.genes.addAll(genesToMiRNA.keySet());
        HashSet<String> outerLoopGenes;
        HashSet<String> innerLoopGenes;
        if(selectedGenes != null){
            outerLoopGenes = selectedGenes;
        }
        else{
            outerLoopGenes = genes;
        }

        //if one gene is selected, consider all other genes as targets
        if(selectedGenes != null && restricted){
            innerLoopGenes = selectedGenes;
        }
        else{
            innerLoopGenes = genes;
        }

        for(String geneA : outerLoopGenes) {

            if(genesToMiRNA.get(geneA) == null){
                logger.warn("No miRNA interactions were found for " + selectedGenes);
                continue;
            }
            for (String geneB : innerLoopGenes) {
                if(geneA.equals(geneB)) continue;

                Object[] sharedMirs = Arrays.stream(genesToMiRNA.get(geneA))
                        .filter(x -> Arrays.stream(genesToMiRNA.get(geneB)).anyMatch(y -> y.equals(x)))
                        .toArray();

                for(Object sharedMiR : sharedMirs){
                    String miR = (String) sharedMiR;
                    this.miRNAs.add(miR);
                    Triplet triple = new Triplet(geneA, geneB, miR);
                    this.triplets.add(triple);
                }
            }
        }
    }

    public void readFileInSetFormat(File fileGenesMiRNA) {
        this.readFileInSetFormat(fileGenesMiRNA, null, false);
    }

    public void readFileInSetFormat(File fileGenesMiRNA, HashSet<String> selectedGenes){
        this.readFileInSetFormat(fileGenesMiRNA, selectedGenes, false);
    }

    /**
     * Reads file with gene - miRNA interactions in the triple format. Initializes miRExpr and geneExpr with the
     * names of genes and miRNAs found in the file.
     */
    public void readFileWithTriples(File fileGenesMiRNA, HashSet<String> selectedGenes){
        readTriples(fileGenesMiRNA, "\t", true, selectedGenes);
    }

    public void readFileWithTriples(File fileGenesMiRNA){
        readTriples(fileGenesMiRNA, "\t", true, null);
    }
}
