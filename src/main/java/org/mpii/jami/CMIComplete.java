package org.mpii.jami;

import java.util.*;

/**
 * Created by fuksova on 11/27/15.
 * For gene pair and one miRNA expression data, it computes CMI and evaluates its p-value via randomization.
 * There are several possibilities how to compute CMI.
 */
public class CMIComplete {
    public static ArrayList<ArrayList<Integer>> randomizedIndices; //used for randomizing data in third array of origData
    public static Random rand;
    ArrayList<double[]> origData;  //Input data for that CMI is computed
    double cmi;  //Obtained cmi value
    double pValue;  //obtained p-value of CMI
    double [] randCMI;
    int numPerm;
    private int maxDeep; //Has an effect only in computeCMIandPValueBetterPartitioning(), sets maxDeep of iterative partitioning grid split


    /**
     * Initializes data for which CMI is computed
     * @param numPerm Number of permutations for p-value computation
     * @param data input data
     */
    public CMIComplete(int numPerm,ArrayList<double[]> data) {
        this.numPerm=numPerm;
        if (randomizedIndices == null||randomizedIndices.size()!=numPerm||randomizedIndices.get(0).size()!=data.get(0).length) {
            initRandomized(data.get(0).length,numPerm);
        }
        if(rand==null){
            rand=new Random();
        }

        this.origData=data;
        maxDeep=Integer.MAX_VALUE;

    }

    /**
     * Initializes random permutations templates
     * @param dataSize length of data in one dimension
     * @param numRandomizations number of randomizations
     */
    public static void initRandomized(int dataSize,int numRandomizations){
        long start=System.currentTimeMillis();
        randomizedIndices=new ArrayList<>(numRandomizations);
        ArrayList<Integer> basic=new ArrayList<>(dataSize);
        for (int i = 0; i < dataSize; i++) {
            basic.add(i);
        }
        for (int i = 0; i < numRandomizations; i++) {
            ArrayList<Integer> newList=new ArrayList<>(basic);
            Collections.shuffle(newList);
            randomizedIndices.add(newList);
        }
        long end=System.currentTimeMillis();
        System.out.println("rand time: "+(end-start));

    }

    /**
     * Old original version of iterative partitioning implemented such that it is as similar to Cupid as possible
     */
    public void computeCMIandPValue(){
        IterativePartitioning ip=new IterativePartitioning(origData);
        randCMI =new double[numPerm];
        cmi=ip.naivePartitioning();
        double[] toBeRandomized=origData.get(2);
        for (int i = 0; i < numPerm; i++) {
            double [] randomizedData=new double[toBeRandomized.length];
            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            for (int j = 0; j < toBeRandomized.length; j++) {
                randomizedData[j]=toBeRandomized[currentPermutation.get(j)];
            }
            ArrayList<double[]> forParitioning=new ArrayList<>(3);
            forParitioning.add(origData.get(0));
            forParitioning.add(origData.get(1));
            forParitioning.add(randomizedData);
            IterativePartitioning ipRand=new IterativePartitioning(forParitioning);
            randCMI[i]=ipRand.naivePartitioning();

        }
        double sum=0;
        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]){
                sum=sum+1.0;
            }
        }
        sum=Math.max(sum,1.0);
        pValue=sum/numPerm;

    }


    /**
     * This function should be used for CMI computation by IterativePartitioning and p-value estimation by randomization
     * of the third data array. This implementation uses the most efficient version of IterativePartitioning.
     * However, drawbacks of this method should be considered. For instance, if all values in the
     * second data set are zeros, this methods outputs nonzero CMI value which can be moreover evaluated as significant.
     */
    public void computeCMIandPValueBetterPartitioning(){
        IterativePartitioning ip=new IterativePartitioning(origData);
        ip.maxDeep=maxDeep;
        randCMI =new double[numPerm];
        cmi=ip.iterativePartitioningBetter();
        double[] toBeRandomized=origData.get(2);
         for (int i = 0; i < numPerm; i++) {

            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            Integer[] randomizedSorted=new Integer[toBeRandomized.length];
            Integer []randomizedInverse=new Integer[toBeRandomized.length];
            for (int j = 0; j < toBeRandomized.length; j++) {

                Integer value = ip.inverseSortedIndices.get(2)[currentPermutation.get(j)];
                randomizedInverse[j]= value;
                randomizedSorted[value]=j;
            }

            ArrayList<Integer[]> sorted=new ArrayList<>();
            sorted.add(ip.sortedIndices.get(0));
            sorted.add(ip.sortedIndices.get(1));
            sorted.add(randomizedSorted);
            ArrayList<Integer[]> inverese=new ArrayList<>();
            inverese.add(ip.inverseSortedIndices.get(0));
            inverese.add(ip.inverseSortedIndices.get(1));
            inverese.add(randomizedInverse);
            IterativePartitioning ipRand=new IterativePartitioning(sorted,inverese);
            ipRand.maxDeep=maxDeep;
            randCMI[i]=ipRand.iterativePartitioningBetter();
//             double control=ipRand.computeCMIInGrid(ipRand.treeRoot);
//             if(Math.abs(randCMI[i]-control)>0.00001){
//                 System.out.println("CMI differs");
//             }
        }
        double sum=0;
        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]){
                sum=sum+1.0;
            }
        }
        sum=Math.max(sum,1.0);
        pValue=sum/numPerm;


    }


    /**
     * This function uses CMI computation on uniform grid
     * @param numberOfBins number of bins of uniform grid in every dimension
     */
    public void computeCMIandPVUniformGrid(int numberOfBins){
        CMIUniform cmiUniform=new CMIUniform(origData);
        randCMI =new double[numPerm];
        cmi=cmiUniform.computeCMI(numberOfBins);
        double[] toBeRandomized=origData.get(2);

        for (int i = 0; i < numPerm; i++) {
            double [] randomized=new double[toBeRandomized.length];
            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            for (int j = 0; j < toBeRandomized.length; j++) {
                randomized[j]=toBeRandomized[currentPermutation.get(j)];
            }
            ArrayList<double[]> randData=new ArrayList<>();
            randData.add(origData.get(0));
            randData.add(origData.get(1));
            randData.add(randomized);
            CMIUniform cmiRand=new CMIUniform(randData);
            randCMI[i]=cmiRand.computeCMI(numberOfBins);
        }
        double sum=0;
        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]){
                sum=sum+1.0;
            }
        }
        sum=Math.max(sum,1.0);
        pValue=sum/numPerm;

    }

    /**
     * This methods combines idea from IterativePartitioning for grid construction with CMI computation onuniform grid.
     * That is, it uses point coordinates orders instead of original double point coordinates. So, instead of bins of
     * equal width, is uses bins with approximately equal point number.
     * @param binNumber number of bins in every dimension
     */
    public void cmiAndPValuePseudoUniform(int binNumber){
        IterativePartitioning ip=new IterativePartitioning(origData);
        randCMI =new double[numPerm];
        cmi=ip.cmiInUniformGrid(binNumber);
        double[] toBeRandomized=origData.get(2);
        for (int i = 0; i < numPerm; i++) {
            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            Integer[] randomizedSorted=new Integer[toBeRandomized.length];
            Integer[] randomizedInverse=new Integer[toBeRandomized.length];
            for (int j = 0; j < toBeRandomized.length; j++) {
                //    randomizedData[j]=toBeRandomized[currentPermutation.get(j)];
                Integer value = ip.inverseSortedIndices.get(2)[currentPermutation.get(j)];
                randomizedInverse[j]= value;
                randomizedSorted[value]=j;
            }

            ArrayList<Integer[]> sorted=new ArrayList<>();
            sorted.add(ip.sortedIndices.get(0));
            sorted.add(ip.sortedIndices.get(1));
            sorted.add(randomizedSorted);
            ArrayList<Integer[]> inverese=new ArrayList<>();
            inverese.add(ip.inverseSortedIndices.get(0));
            inverese.add(ip.inverseSortedIndices.get(1));
            inverese.add(randomizedInverse);
            IterativePartitioning ipRand=new IterativePartitioning(sorted,inverese);

            randCMI[i]=ipRand.cmiInUniformGrid(binNumber);
        }
        double sum=0;
        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]){
                sum=sum+1.0;
            }
        }
        sum=Math.max(sum,1.0);
        pValue=sum/numPerm;

    }


    /**
     * This methods uses original Iterative Partitioning method for CMI computation. However, for p-value evaluation,
     * CMI grid obtained on first randomization of data is used for all other randomizations.
     */
    public void cmiAndPValueFixedGrid(){
        IterativePartitioning ip=new IterativePartitioning(origData);
        randCMI =new double[numPerm];

        cmi=ip.iterativePartitioningBetter();
        double[] toBeRandomized=origData.get(2);
        TreeNode grid = null;
        LinkedList<Cube> fixedGrid=null;
       // LinkedList<Cube> fixedGrid=ip.usedCubes;
        for (int i = 0; i < numPerm; i++) {
            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            Integer[] randomizedSorted=new Integer[toBeRandomized.length];
            Integer []randomizedInverse=new Integer[toBeRandomized.length];
            for (int j = 0; j < toBeRandomized.length; j++) {

                Integer value = ip.inverseSortedIndices.get(2)[currentPermutation.get(j)];
                randomizedInverse[j]= value;
                randomizedSorted[value]=j;
            }

            ArrayList<Integer[]> sorted=new ArrayList<>();
            sorted.add(ip.sortedIndices.get(0));
            sorted.add(ip.sortedIndices.get(1));
            sorted.add(randomizedSorted);
            ArrayList<Integer[]> inverese=new ArrayList<>();
            inverese.add(ip.inverseSortedIndices.get(0));
            inverese.add(ip.inverseSortedIndices.get(1));
            inverese.add(randomizedInverse);

            IterativePartitioning ipRand=new IterativePartitioning(sorted,inverese);

            if(grid==null){
                randCMI[i]=ipRand.iterativePartitioningBetter();
                grid=ipRand.treeRoot;
                fixedGrid=ipRand.usedCubes;
            }
            else{
                randCMI[i]=ipRand.computeCMIInGrid(grid,fixedGrid);

            }
        }
        double sum=0;

        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]){
                sum=sum+1.0;
            }
        }
        sum=Math.max(sum,1.0);

        pValue=sum/numPerm;

    }

    /**
     * This mediator status is based on Cupid method. It computes pValue^weight, where weight is based on step 2 score.
     * If this value is smaller than certain threshold given by cutoff, it is set to true, which means that this
     * result should be selected as significant.
     * @param weight weight computed in SelectMiRNACandidates
     * @param cutoff threshold on pValue^weight
     * @return significant or not
     */
    public boolean mediatorStatus(double weight, double cutoff){
        double mediatorScore=Math.pow(pValue,weight);
        if(mediatorScore<cutoff){
            return true;
        }
        else{
            return false;
        }
    }

    /**
     * Has an effect only in computeCMIandPValueBetterPartitioning(), sets maxDeep of iterative partitioning grid split
     * @param maxDeep
     */
    public void setMaxDeep(int maxDeep) {
        this.maxDeep = maxDeep;
    }

    public int getMaxDeep() {
        return maxDeep;
    }
}
