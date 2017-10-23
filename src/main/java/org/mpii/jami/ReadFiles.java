package org.mpii.jami;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.IllegalFormatException;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 11/19/15.
 */
public class ReadFiles {

    public static int geneInd;
    public static int miRNAInd;
    public static int valueInd;

    /**
     * Extract gene/miRNA identifiers from a file
     * @param file name of the file from which we want extract the names
     * @param separator How the columns in the file are separated
     * @param indices The array should contain indices of columns where the names are stored
     * @return extracted names
     */
    public static ArrayList<String> extractNamesFromColumn(File file, String separator, boolean skipFirst, int[] indices){
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        ArrayList<String> names=new ArrayList<>();
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

                for (int i = 0; i < indices.length; i++) {
                    names.add(entries[indices[i]]);
                }

                line=br.readLine();
            }
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return names;
    }

    /**
     * Extract names (gene/miRNA identifiers)  from a file.
     * @param file Name of file
     * @param separator separator of columns
     * @param lineNumber number of line where names are stored
     * @return list of names
     */
    public static ArrayList<String> extractNamesFromLine(File file,String separator,int lineNumber){
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        ArrayList<String> names=new ArrayList<>();
        try {
            fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line=br.readLine();
            int counter=0;
            while(counter<lineNumber){
                line=br.readLine();
                counter++;
            }
                String[] entries = pattern.split(line);
                for (int i = 0; i < entries.length; i++) {
                    names.add(entries[i]);
                }


            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return names;
    }

    /**
     * Simplified version of extracting names assumes that all names are stored in the first column of the file
     * @param file name of the file from which we want extract the names
     * @param separator How are the columns in the file separated
     * @return extracted names
     */
    public static ArrayList<String> extractNamesFromColumn(File file, String separator, boolean skipFirst){
        int[] indices=new int[1];
        return extractNamesFromColumn(file, separator, skipFirst, indices);
    }

    /**
     * File for counting lines in a file. It returns the values as long. Needs to be checked once again. Maybe it does not
     * output correct value.
     * @param file name of file
     * @return number of lines in the file
     */
    public static long getLineNumber(File file){
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(file));
            lnr.skip(Long.MAX_VALUE);
            long lineNumber=(lnr.getLineNumber() + 1); //Maybe without +1?
            lnr.close();
            return lineNumber;
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return -1;

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
    public static HashMap<String,String[]> geneToMiRNA(File file,String separatorMain,String separatorMiRNA,boolean skipFirst){
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
     * Reads file with triple of two genes and one miRNA names whose interactions will be tested.
     * Every row contains one triple to test separated by "separator".
     * @param file Name of file
     * @param separator Column separator
     * @param skipFirst Skip first row
     * @return ArrayList with all triple from the file
     */
    public static ArrayList<String[]> readTriples(File file,String separator,boolean skipFirst){
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

        return result;
    }






}
