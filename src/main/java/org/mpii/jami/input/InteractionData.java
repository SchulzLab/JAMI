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

    public ArrayList<Triplet> getTriplets() {
        return triplets;
    }

    private ArrayList<Triplet> triplets = new ArrayList<>();

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
    private void readTriples(File file, String separator, boolean skipFirst, String selectedGene){
        Pattern pattern = Pattern.compile(separator);

        try {
            BufferedReader br = readFile(file);
            String line;
            if(skipFirst){
                br.readLine();
            }
            line=br.readLine();
            while (line!=null) {
                String[] entries = pattern.split(line);
                if((selectedGene != null && entries[0].equals(selectedGene)) || selectedGene == null) {
                    Triplet entry = new Triplet(entries);
                    if (!this.triplets.contains(entry)) this.triplets.add(entry);

                    this.genes.add(entries[0]);
                    this.genes.add(entries[1]);
                    this.miRNAs.add(entries[2]);
                }
                if(selectedGene == null) {
                    Triplet reverseEntry = new Triplet(entries[1], entries[0], entries[2]);
                    if (!this.triplets.contains(reverseEntry)) this.triplets.add(reverseEntry);
                }

                line=br.readLine();
            }
            br.close();

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
    public void readFileInSetFormat(File fileGenesMiRNA, String selectedGene){
        HashMap<String, String[]> genesToMiRNA = geneToMiRNA(fileGenesMiRNA, "\t", ",", true);

        this.genes.addAll(genesToMiRNA.keySet());
        HashSet<String> outerLoopGenes;
        if(selectedGene != null){
            outerLoopGenes = new HashSet<>();
            outerLoopGenes.add(selectedGene);
        }
        else{
            outerLoopGenes = genes;
        }

        for(String geneA : outerLoopGenes) {

            for (String geneB : genes) {
                if(geneA.equals(geneB)) continue;

                HashSet<String> miRNAsA =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneA)));

                HashSet<String> miRNAsB =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneB)));

                miRNAsA.retainAll(miRNAsB); //miRNAsA contains only shared miRNAs
                this.miRNAs.addAll(miRNAsA);

                for (String miR : miRNAsA) {
                    Triplet triple = new Triplet(geneA, geneB, miR);
                    if(!this.triplets.contains(triple)) this.triplets.add(triple);
                }
            }
        }
    }

    /**
     * Reads file with gene - miRNA interactions in the triple format. Initializes miRExpr and geneExpr with the
     * names of genes and miRNAs found in the file.
     */
    public void readFileWithTriples(File fileGenesMiRNA, String selectedGene){
        readTriples(fileGenesMiRNA, "\t", true, selectedGene);
    }

    public void filterByGene(String selectedGene) {
        ArrayList<Triplet> retainedTriplets = new ArrayList<>();
        for(Triplet triplet : triplets){
            if(triplet.getGeneOne().equals(selectedGene))
                retainedTriplets.add(triplet);
        }
        this.triplets = retainedTriplets;
    }
}
