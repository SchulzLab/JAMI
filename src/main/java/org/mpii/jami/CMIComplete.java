package org.mpii.jami;

import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

/**
 * Created by fuksova on 11/27/15.
 * For gene pair and one miRNA expression data, it computes CMI and evaluates its p-value via randomization.
 * There are several possibilities how to compute CMI.
 */
public class CMIComplete{
    public static ArrayList<ArrayList<Integer>> randomizedIndices; //used for randomizing data in third array of origData
    ArrayList<double[]> origData;  //Input data for that CMI is computed
    double cmi;  //Obtained cmi value
    double pValue;  //obtained p-value of CMI
    int numPerm; //number of permutations for p-value computation
    private int maxDeep; //maximum depth of iterative partitioning grid split
    ForkJoinPool commonPool;


    /**
     * Initializes data for which CMI is computed
     * @param numPerm Number of permutations for p-value computation
     * @param data input data
     */
    public CMIComplete(int numPerm,ArrayList<double[]> data, ForkJoinPool fjpool) {
        this.numPerm=numPerm;
        this.origData=data;
        maxDeep=Integer.MAX_VALUE;
        commonPool = fjpool;
    }

    /**
     * Initializes random permutations templates
     * @param dataSize length of data in one dimension
     * @param numRandomizations number of randomizations
     */
    public static void initRandomized(int dataSize,int numRandomizations){
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
    }

    /**
     * Old original version of iterative partitioning implemented such that it is as similar to Cupid as possible
     */
    public void computeAsCUPID(){
        IterativePartitioning ip=new IterativePartitioning(origData);
        double[] randCMI =new double[numPerm];
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
    public void computeIterativePartitioning(){
        IterativePartitioning ip=new IterativePartitioning(origData);
        ip.maxDeep=maxDeep;
        cmi=ip.iterativePartitioningBetter();

        RecursiveTask<Integer> task = new CMIRecursiveTask(numPerm, cmi, ip, origData.get(2).length);
        int sum = commonPool.invoke(task);
        sum=Math.max(sum, 1);
        pValue=(double) sum/numPerm;
    }

    /**
     * This function uses CMI computation on uniform grid
     * @param numberOfBins number of bins of uniform grid in every dimension
     */
    public void computeUniformGrid(int numberOfBins){
        CMIUniform cmiUniform=new CMIUniform(origData);
        cmi=cmiUniform.computeCMI(numberOfBins);

        RecursiveTask<Integer> task = new CMIUniformGridRecursiveTask(numPerm, cmi, origData, origData.get(2).length, numberOfBins);
        int sum = commonPool.invoke(task);
        sum=Math.max(sum, 1);
        pValue=(double) sum/numPerm;
    }

    /**
     * This methods combines idea from IterativePartitioning for grid construction with CMI computation on uniform grid.
     * That is, it uses point coordinates orders instead of original double point coordinates. So, instead of bins of
     * equal width, is uses bins with approximately equal point number.
     * @param numberOfBins number of bins in every dimension
     */
    public void computePseudoUniformGrid(int numberOfBins){
        IterativePartitioning ip=new IterativePartitioning(origData);

        cmi=ip.cmiInUniformGrid(numberOfBins);

        RecursiveTask<Integer> task = new CMIUniformGridRecursiveTask(numPerm, cmi, origData, origData.get(2).length, numberOfBins);
        int sum = commonPool.invoke(task);
        sum=Math.max(sum, 1);
        pValue=(double) sum/numPerm;
    }
}
