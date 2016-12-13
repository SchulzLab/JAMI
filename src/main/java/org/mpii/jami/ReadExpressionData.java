package org.mpii.jami;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 12/1/15.
 * Methods for reading expression data files and storing those data in suitable format.
 * There is one specific integer value for every gene/miRNA. This is useful for finding corresponding expression
 * data in the field expressionData.
 */
public class ReadExpressionData {

    private ArrayList<String> names; //Names stored at corresponding indices
    private ArrayList<double[]> expressionData; //expressionData.get(i) returns expression data corresponding to the i-th name
    private HashMap<String,Integer> nameToData; //Returns index corresponding to given name
    boolean genesAreRows;
    int numberOfSamples;
    boolean skipFirst;
    double multiplier;
    boolean modifyName;
    String separator;


    /**
     * In constructor only fields for storing names are initialized. Those fields can be changed during
     * methods for file reading if there is no data for given name.
     * @param names List of names for which the expression data should be extracted
     */
    public ReadExpressionData(ArrayList<String> names){
        this.names=names;
        nameToData=new HashMap<>();
        for (int i = 0; i < names.size(); i++) {
            nameToData.put(names.get(i),i);
        }
         genesAreRows=true;
         numberOfSamples =-1;
         skipFirst=true;
         multiplier=1;
         modifyName=false;
        separator="\t";
    }


    /**
     * Reads the file and stores data
     * @param fileName name of the file
     */
    public void readFile(String fileName){
        if(genesAreRows){
            readGenesAreRows(fileName);
        }
        else{
            if(numberOfSamples >=1){
                readGenesAreColumns(fileName, numberOfSamples);
            }
            else{
                readGenesAreColumns(fileName);
            }
        }
    }


    /**
     * Reads file in such an input format where one row contain expression data for one certain gene/miRNA and one column
     * contains all gene/miRNA expression data for one sample.
     * @param fileName Name of file with gene expression data
      */
    private void readGenesAreRows(String fileName){ //Add boolean modifyName
        expressionData = new ArrayList<>(names.size());
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        try {
            fr = new FileReader(fileName);
            BufferedReader br = new BufferedReader(fr);

            if (names.isEmpty()) {
                throw new IllegalArgumentException("Array names can not be empty");
            }
            for (int i = 0; i < names.size(); i++) {
                expressionData.add(null);
            }
            HashSet<String> notUsedNames = new HashSet<>();
            notUsedNames.addAll(names);
            if(skipFirst){
                br.readLine(); //Skips first line
            }
            String line = br.readLine();
            while (line != null&(!notUsedNames.isEmpty())) {
                String[] entries = pattern.split(line);
                String extractedName;
                if(modifyName){
                    extractedName= modifyName(entries[0]);
                }
                else{
                    extractedName=entries[0];
                }
                if (notUsedNames.contains(extractedName)) {
                    int countOfSamples = entries.length - 1;
                    double[] doubleEntries = new double[countOfSamples];
                    for (int i = 1; i <= countOfSamples; i++) {
                        double entry = Double.parseDouble(entries[i]);
                        entry=entry*multiplier;
                        doubleEntries[i-1]=entry;
                    }
                    notUsedNames.remove(extractedName);
                    expressionData.set(nameToData.get(extractedName), doubleEntries);

                }
                line = br.readLine();
            }
            ArrayList<String> newNames=new ArrayList<>();
            ArrayList<double[]> newExpressionData=new ArrayList<>();
            nameToData.clear();
            for (int i = 0; i < names.size(); i++) {
                String name=names.get(i);
                if(!notUsedNames.contains(name)){
                    newNames.add(name);
                    newExpressionData.add(expressionData.get(i));
                    nameToData.put(name,newNames.size()-1);
                }
            }
            names=newNames;
            expressionData=newExpressionData;

            br.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

    /**
     * From given string remove quotes and cut all characters behind the first occurrence of a dot
     * @param identifier
     * @return
     */
    public String modifyName(String identifier){
        String name=identifier.replaceAll("\"","");
        int endIndex = name.indexOf(".");
        String name2;
        if(endIndex>0){
            name2=name.substring(0, endIndex);
        }
        else{
            name2=name;
        }

        return name2;
    }


    /**
     * Reads file in such an input format where one column contain expression data for one certain gene/miRNA and one
     * row contains all gene/miRNA expression data for one sample.
     * @param fileName Name of file with gene expression data
     * @param countOfSamples How many rows with numerical data does the file contain
     *                        numerical values, otherwise an error occurs.
     */

    private void readGenesAreColumns(String fileName,int countOfSamples){
        expressionData = new ArrayList<>(names.size());
        Pattern pattern = Pattern.compile(separator);
        FileReader fr;
        int temp=0;
        if(skipFirst){
            temp=1;
        }
        try {
            fr = new FileReader(fileName);
            BufferedReader br = new BufferedReader(fr);

            if (names.isEmpty()) {
                throw new IllegalArgumentException("Array names can not be empty");
            }
            for (int i = 0; i < names.size(); i++) {
                expressionData.add(new double[countOfSamples]);
            }
            ArrayList<Integer> columnIndices=new ArrayList<>(names.size());
            for (int i = 0; i < names.size(); i++) {
                columnIndices.add(-1);
            }
            String line = br.readLine();
            String[] entries = pattern.split(line);
            for (int i = 0; i < entries.length; i++) {//It should be enough for the columnIndices to have length of names.size
                String extractedName;
                if(modifyName){
                    extractedName= modifyName(entries[i]);
                }
                else{
                    extractedName=entries[i];
                }
                Integer index=nameToData.get(extractedName);
                if(index!=null){
                    columnIndices.set(index,i+temp);  //the name with index in the list corresponds (i+1)-th column in file (first column contains sample ids)
                }
            }
            line=br.readLine();
            int sampleCounter=0;
            while (line != null&&sampleCounter<countOfSamples) {
                entries = pattern.split(line);
                for (int i = 0; i < columnIndices.size(); i++) {
                    int index=columnIndices.get(i);  //gives index in entries of value corresponding to i-th name
                    if(index!=-1){
                        Double value=Double.parseDouble(entries[index]);
                        value=value*multiplier;
                        expressionData.get(i)[sampleCounter]=value;
                    }
//                    else if(sampleCounter==0){
//                        System.out.println("Gene with name "+names.get(i)+" "+"is missing");
//                    }


                }
                line = br.readLine();
                sampleCounter++;
            }//create new names, new expression data new hash map
            ArrayList<String> newNames=new ArrayList<>();
            ArrayList<double[]> newExpressionData=new ArrayList<>();
            nameToData.clear();
            for (int i = 0; i < columnIndices.size(); i++) {
                Integer index = columnIndices.get(i);
                if(index !=-1){
                    newNames.add(names.get(i));
                    newExpressionData.add(expressionData.get(i));
                    nameToData.put(names.get(i),newNames.size()-1);
                }
            }

            names=newNames;
            expressionData=newExpressionData;

            br.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }

    }

    private void readGenesAreColumns(String fileName){
        int lineNumber=(int)ReadFiles.getLineNumber(fileName);//seems to count one more
        if(skipFirst) lineNumber--;
        readGenesAreColumns(fileName,lineNumber);

    }

    /**
     * for debugging purposes: Checks whether name<->index fields are correct
     */
    public void checkIndices(){
        for (int i = 0; i < names.size(); i++) {
            String name=names.get(i);
            int index=nameToData.get(name);
            if(index!=i){
                System.out.println("Error in indices. Name index "+i);
            }
        }
    }

    /**
     * Name of gene/miRNA at index i means that expression data to this name are stored at the i-th index in expressionData
     * @return list with gene/miRNA names
     */
    public ArrayList<String> getIntegersToNames() {
        return names;
    }

    public ArrayList<double[]> getExpressionData() {
        return expressionData;
    }

    public HashMap<String, Integer> getNameToData() {
        return nameToData;
    }

    /**
     * Expression data to gene/miRNA name
     * @param id gene/miRNA name
     * @return expression data for given gene/miRNA name
     */
    public double[] getExpressionDataToID(String id){
        Integer index=nameToData.get(id);
        if(index!=null){
            return expressionData.get(index);
        }
        else{
            return null;
        }
    }

}
