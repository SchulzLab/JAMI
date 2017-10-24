package org.mpii.jami;

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

    private HashSet<String> genes = new HashSet<>();
    private HashSet<String> miRNAs = new HashSet<>();

    public HashSet<String> getGenes() {
        return genes;
    }

    public HashSet<String> getMiRNAs() {
        return miRNAs;
    }

    public ArrayList<String[]> getTriplets() {
        return triplets;
    }

    private ArrayList<String[]> triplets;

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
    public void readTriples(File file, String separator, boolean skipFirst){
        Pattern pattern = Pattern.compile(separator);
        ArrayList<String[]> result=new ArrayList<>();
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
                result.add(entries);
                line=br.readLine();
            }
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        this.triplets = result;
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
                    br.readLine();
                    skipFirst = false;
                    continue;
                }
                String[] entries = pattern.split(line);

                if(entries.length!=2){  //maybe add possibility that there is no miRNA
                    throw new IllegalArgumentException("Wrong line format in file "+file);
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
    protected void readFileInSetFormat(File fileGenesMiRNA){
        HashMap<String, String[]> genesToMiRNA = geneToMiRNA(fileGenesMiRNA, "\t", ",", false);

        this.genes.addAll(genesToMiRNA.keySet());

        for (String[] miRNANames : genesToMiRNA.values()) {
            for (int i = 0; i < miRNANames.length; i++) {
                miRNAs.add(miRNANames[i]);
            }
        }

        for(String geneA : genes) {
            for (String geneB : genes) {
                if (geneA == geneB) continue; //do not test same gene

                HashSet<String> miRNAsA =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneA)));

                HashSet<String> miRNAsB =
                        new HashSet<>(Arrays.asList(genesToMiRNA.get(geneB)));

                miRNAsA.retainAll(miRNAsB); //miRNAsA contains only shared miRNAs

                for (String miR : miRNAsA) {
                    String[] triple = {geneA, geneB, miR};
                    this.triplets.add(triple);
                }
            }
        }
    }

    /**
     * Reads file with gene - miRNA interactions in the triple format. Initializes miRExpr and geneExpr with the
     * names of genes and miRNAs found in the file.
     */
    protected void readFileWithTriples(File fileGenesMiRNA){
        readTriples(fileGenesMiRNA, "\t", false);
        for (String[] stringTriple : this.triplets) {
            this.genes.add(stringTriple[0]);
            this.genes.add(stringTriple[1]);
            this.miRNAs.add(stringTriple[2]);
        }
    }
}
