package org.mpii.jami.scores;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by fuksova on 11/20/15.
 */
public class ScoreTable {
    //For merging scores of two genes
    //Should be replaced by Score Structure, resp. SparseScoreMatrix
    private ArrayList<String> miRNAnames;
    double[][] data;

    /**
     * Creates table with score values for two genes and a set of miRNA
     * @param gene1Data map: miRNA name -> score with gene 1
     * @param gene2Data map: miRNA name -> score with gene 2
     */
    public ScoreTable(HashMap<String,Double> gene1Data,HashMap<String,Double> gene2Data){
        ArrayList<Double> firstGene=new ArrayList<>();
        ArrayList<Double> secondGene=new ArrayList<>();
        miRNAnames=new ArrayList<>();
        for (Map.Entry<String, Double> entry : gene1Data.entrySet()) {
            String key = entry.getKey();
            Double value=entry.getValue();
            Double value2 = gene2Data.remove(key);
            if(value2==null){
                value2=0.0;
            }
            miRNAnames.add(key);
            firstGene.add(value);
            secondGene.add(value2);
        }
        for (Map.Entry<String, Double> entry : gene2Data.entrySet()) {
            String key = entry.getKey();
            Double value=entry.getValue();
            Double value2 = gene1Data.remove(key);
            if(value2==null){
                value2=0.0;
            }
            miRNAnames.add(key);
            firstGene.add(value2);
            secondGene.add(value);
        }

        data=new double[firstGene.size()][2];
        for (int i = 0; i < firstGene.size(); i++) {
            data[i][0]=firstGene.get(i);
            data[i][1]=secondGene.get(i);
        }
    }

    /**
     * Creates table with score values for two genes and a set of miRNA whose names are given in names in range start
     * to end
     * @param gene1Data map: miRNA name -> score with gene 1
     * @param gene2Data map: miRNA name -> score with gene 2
     * @param names list with names of miRNAs
     * @param start index in names where to start inclusive
     * @param end index in names where to end exclusive
     */
    public ScoreTable(HashMap<String,Double> gene1Data,HashMap<String,Double> gene2Data,ArrayList<String> names,int start,int end){
        miRNAnames=new ArrayList<>(names.size());
        data=new double[names.size()][2];
        for (int i = start; i < end; i++) {
            String name=names.get(i);
            miRNAnames.add(name);
            Double value1 = gene1Data.get(name);
            Double value2 = gene2Data.get(name);
            if(value2==null){
                value2=0.0;
            }
            if(value1==null){
                value1=0.0;
            }
            data[i-start][0]=value1;
            data[i-start][1]=value2;
        }

    }

    public double[][] getData() {
        return data;
    }

    public ArrayList<String> getMiRNAnames() {
        return miRNAnames;
    }
}
