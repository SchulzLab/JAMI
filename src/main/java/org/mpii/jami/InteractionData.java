package org.mpii.jami;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.mpii.jami.model.Triplet;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

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
    private void readTriples(File file, String separator, boolean skipFirst){
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;

        try {
            fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line;
            if(skipFirst){
                br.readLine();
            }
            line=br.readLine();
            while (line!=null) {
                String[] entries = pattern.split(line);
                this.triplets.add(new Triplet(entries));
                this.genes.add(entries[0]);
                this.genes.add(entries[1]);
                this.miRNAs.add(entries[2]);
                line=br.readLine();
            }
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
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
        FileReader fr;
        HashMap<String,String[]> genesWithMiRNA=new HashMap<>();
        try {
            fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
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
    public void readFileInSetFormat(File fileGenesMiRNA){
        HashMap<String, String[]> genesToMiRNA = geneToMiRNA(fileGenesMiRNA, "\t", ",", true);

        this.genes.addAll(genesToMiRNA.keySet());
        ArrayList<String> genes_processed = new ArrayList<>();

        for(String geneA : genes) {
            genes_processed.add(geneA);
            ArrayList<String> genes_to_test = new ArrayList<>(genes);
            genes_to_test.removeAll(genes_processed);

            for (String geneB : genes_to_test) {

                HashSet<String> miRNAsA =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneA)));

                HashSet<String> miRNAsB =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneB)));

                miRNAsA.retainAll(miRNAsB); //miRNAsA contains only shared miRNAs
                this.miRNAs.addAll(miRNAsA);

                for (String miR : miRNAsA) {
                    Triplet triple = new Triplet(geneA, geneB, miR);
                    this.triplets.add(triple);
                }
            }
        }
    }

    /**
     * Reads file with gene - miRNA interactions in the triple format. Initializes miRExpr and geneExpr with the
     * names of genes and miRNAs found in the file.
     */
    public void readFileWithTriples(File fileGenesMiRNA){
        readTriples(fileGenesMiRNA, "\t", true);
    }

    public void filterByGene(String selectedGene) {
        ArrayList<Triplet> retainedTriplets = new ArrayList<>();
        for(Triplet triplet : triplets){
            if(triplet.getGeneOne().equals(selectedGene))
                retainedTriplets.add(triplet);
            if(triplet.getGeneTwo().equals(selectedGene)){
                retainedTriplets.add(new Triplet(triplet.getGeneTwo(), triplet.getGeneOne(), triplet.getMiRNA()));
            }
        }
        this.triplets = retainedTriplets;
    }
}
