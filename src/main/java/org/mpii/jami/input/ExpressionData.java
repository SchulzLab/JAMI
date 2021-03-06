package org.mpii.jami.input;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.core.util.FileUtils;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by fuksova on 12/1/15.
 * Methods for reading expression data files and storing those data in suitable format.
 * There is one specific integer value for every gene/miRNA. This is useful for finding corresponding expression
 * data in the field expressionData.
 */
public class ExpressionData {

    private ArrayList<String> names; //Names stored at corresponding indices
    private ArrayList<List<Double>> expressionData; //expressionData.get(i) returns expression data corresponding to the i-th name
    private HashMap<String,Integer> nameToData; //Returns index corresponding to given name
    private int numberOfSamples = -1;

    public int getNumberOfSamples() {
        return numberOfSamples;
    }

    /**
     * In constructor only fields for storing names are initialized. Those fields can be changed during
     * methods for file reading if there is no data for given name.
     * @param names List of names for which the expression data should be extracted
     */
    public ExpressionData(ArrayList<String> names){
        this.names = names;
        nameToData=new HashMap<>();
        expressionData = new ArrayList<>();
    }

    public ExpressionData(){
        nameToData=new HashMap<>();
        expressionData = new ArrayList<>();
    }

    /**
     * Reads the file and stores data
     * @param file name of the file
     */
    public void readFile(File file, boolean hasHeader){
        try {
            int entryCount = 0;

            CSVFormat format;
            if(hasHeader) format = CSVFormat.TDF.withFirstRecordAsHeader();
            else format = CSVFormat.TDF;

            CSVParser parser;

            if(FileUtils.getFileExtension(file).equals("gz"))
            {
                FileInputStream fis = new FileInputStream(file);
                GZIPInputStream gis = new GZIPInputStream(fis);
                InputStreamReader isr = new InputStreamReader(gis);
                parser = CSVParser.parse(isr, format);
            }
            else parser = CSVParser.parse(file, Charset.defaultCharset(), format);

            for (CSVRecord record : parser)
            {
                if(this.numberOfSamples == -1)
                    this.numberOfSamples = record.size() - 1;
                if(record.size()-1 != this.numberOfSamples || this.numberOfSamples < 1){
                    throw new IOException(" record " + record.getRecordNumber() +
                            " has the wrong number of elements." );
                }
                if(names != null){
                    if(names.contains(record.get(0))) {
                        addEntry(record, entryCount++); //add only selected entries
                    }
                } else{
                    addEntry(record, entryCount++); //adding all entries
                }
            }
        }catch(FileNotFoundException fnfe){
            System.err.println("File" + file.getAbsolutePath() + "not found.");
            System.exit(1);
        }catch(IOException ioe){
            System.err.println("IO error" + ioe.getMessage());
            System.exit(1);
        }
    }

    private void addEntry(CSVRecord record, int entryCount){
        nameToData.put(record.get(0), entryCount);
        ArrayList<Double> entries = new ArrayList<>(record.size()-1);
        for (int i = 1; i < record.size(); i++) {
            entries.add(Double.parseDouble(record.get(i)));
        }
        expressionData.add(Collections.unmodifiableList(entries));
    }

    /**
     * Name of gene/miRNA at index i means that expression data to this name are stored at the i-th index in expressionData
     * @return list with gene/miRNA names
     */

    public ArrayList<List<Double>> getExpressionData() {
        return expressionData;
    }

    public HashMap<String, Integer> getNameToData() {
        return nameToData;
    }


}
