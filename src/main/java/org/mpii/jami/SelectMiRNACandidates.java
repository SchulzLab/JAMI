package org.mpii.jami;

import java.util.ArrayList;

/**
 * Created by fuksova on 11/19/15.
 * From a set of miRNA, select those that should be inspected together with given gene pair. Based on Cupid scores and
 * methods
 */
public class SelectMiRNACandidates {

    private ArrayList<String> names;  //list of candidate miRNA names
    int numberOfPermutations;  //Number of permutations that will be used in evaluation in CMI significance
    private double [][] score; //score obtained from the score file, first index: miRNA, second index: 0: first gene, 1: second gene
    private double [] selectedWeights; //weights computed as geomean of score
    private int[] selectedIndices; //Indices of selected miRNA. Be carefull, this filtering is not sufficient,
                // even miRNA with constant zero expression values can pass
    private int numberSelected; //Number of miRNA that were selected


    /**
     *
     * @param names List with miRNA names
     * @param numberOfPermutations Number of permutations that will be used in evaluation in CMI significance
     * @param score score[i][0] score of the interaction of the i-th miRNA with gene in data.get(0)
     */
    public SelectMiRNACandidates(ArrayList<String> names, int numberOfPermutations, double[][] score){
        this.names=names;
        this.numberOfPermutations=numberOfPermutations;
        this.score=score;

    }

    /**
     * Select candidate miRNA from the whole list. Condition is p^sqrt(S_im*S_jm)>cutoff. Where S_im is the score for
     * i-th gene and m-th miRNA
     * @param cutoff score cutoff
     */
    public void selectCandidates(double cutoff){

        double [] weights=geomean();
        double threshold=Math.log10(cutoff)/Math.log10(1.0/numberOfPermutations);
        ArrayList<Integer> selected=new ArrayList<>();
        for (int i = 0; i < score.length; i++) {
            if(weights[i]>threshold){
                selected.add(i);
            }
        }
        numberSelected = selected.size();
        selectedIndices=new int[numberSelected];
        for (int i = 0; i < numberSelected; i++) {
            selectedIndices[i]=selected.get(i);
        }
        selectedWeights = new double[numberSelected];
        for (int i = 0; i < numberSelected; i++) {
            selectedWeights[i]=weights[selectedIndices[i]];
        }


    }


    /**
     * Geometric mean of an array with doubles
     * @return geometric mean
     */
    public double[] geomean(){
        double result[]=new double[score.length];
        int l=score[0].length;
        double pow=1.0/l;
        for (int i = 0; i < score.length; i++) {
            double temp=1;
            for (int j = 0; j < l; j++) {
                temp*=score[i][j];
            }
            if(l==2) {
                result[i] = Math.sqrt(temp);
            }
            else{
                result[i]=Math.pow(temp, pow);  //maybe change to sqrt for two dimensional data maybe it has a better precision or something
            }
        }
        return result;
    }

    /**
     *
     * @return weights corresponding to the selected miRNA
     */
    public double[] getSelectedWeights() {
        return selectedWeights;
    }

    /**
     *
     * @return score obtained from the score file, first index: miRNA, second index: 0: first gene, 1: second gene
     */
    public double[][] getScore() {
        return score;
    }

    /**
     *
     * @return miRNA names
     */
    public ArrayList<String> getNames() {
        return names;
    }

    /**
     *
     * @return indices of selected miRNA
     */
    public int[] getSelectedIndices() {
        return selectedIndices;
    }

    /**
     *
     * @return number of selected miRNA for this gene pair
     */
    public int getNumberSelected() {
        return numberSelected;
    }
}
