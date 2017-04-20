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
     * @param fileName name of the file from which we want extract the names
     * @param separator How the columns in the file are separated
     * @param indices The array should contain indices of columns where the names are stored
     * @return extracted names
     */
    public static ArrayList<String> extractNamesFromColumn(String fileName, String separator, boolean skipFirst, int[] indices){
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        ArrayList<String> names=new ArrayList<>();
        try {
            fr = new FileReader(fileName);
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
     * @param fileName Name of file
     * @param separator separator of columns
     * @param lineNumber number of line where names are stored
     * @return list of names
     */
    public static ArrayList<String> extractNamesFromLine(String fileName,String separator,int lineNumber){
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        ArrayList<String> names=new ArrayList<>();
        try {
            fr = new FileReader(fileName);
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
     * @param fileName name of the file from which we want extract the names
     * @param separator How are the columns in the file separated
     * @return extracted names
     */
    public static ArrayList<String> extractNamesFromColumn(String fileName, String separator, boolean skipFirst){
        int[] indices=new int[1];
        return extractNamesFromColumn(fileName, separator, skipFirst, indices);
    }

    /**
     * File for counting lines in a file. It returns the values as long. Needs to be checked once again. Maybe it does not
     * output correct value.
     * @param fileName name of file
     * @return number of lines in the file
     */
    public static long getLineNumber(String fileName){
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(new File(fileName)));
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


    public static HashMap<String,String[]> geneToMiRNA(String fileName,String separatorMain,String separatorMiRNA,boolean skipFirst){
        Pattern pattern = Pattern.compile(separatorMain);
        Pattern patternMiRNA=Pattern.compile(separatorMiRNA);
        FileReader fr;
        HashMap<String,String[]> genesWithMiRNA=new HashMap<>();
        try {
            fr = new FileReader(fileName);
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
                    throw new IllegalArgumentException("Wrong line format in file "+fileName);
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


    public static ArrayList<String[]> readTriples(String fileName,String separatorMain,boolean skipFirst){
        Pattern pattern = Pattern.compile(separatorMain);
        ArrayList<String[]> result=new ArrayList<>();
        FileReader fr;

        try {
            fr = new FileReader(fileName);
            BufferedReader br = new BufferedReader(fr);
            String line=br.readLine();
            while (line!=null) {
                if(skipFirst){
                    br.readLine();
                    skipFirst = false;
                    continue;
                }
                String[] entries = pattern.split(line);
                result.add(entries);
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
