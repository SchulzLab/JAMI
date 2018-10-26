package org.mpii.jami;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.mpii.jami.cmi.CMIUniform;
import org.mpii.jami.cmi.Cube;
import org.mpii.jami.cmi.IterativePartitioning;
import org.mpii.jami.helpers.ComputeCMI;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

/**
 * Created by fuksova on 11/27/15.
 * For gene pair and one miRNA expression data, it computes CMI and evaluates its p-value via randomization.
 * There are several possibilities how to compute CMI.
 */
public class CMIComplete{
    public static ArrayList<ArrayList<Integer>> randomizedIndices; //used for randomizing data in third array of origData
    ArrayList<List<Double>> origData;  //Input data for that CMI is computed
    double cmi;  //Obtained cmi value
    double pValue;  //obtained p-value of CMI
    int numPerm; //number of permutations for p-value computation
    private int maxDeep; //maximum depth of iterative partitioning grid split
    private boolean considerZeros;
    private Random random;

    /**
     * Initializes data for which CMI is computed
     * @param numPerm Number of permutations for p-value computation
     * @param data input data
     */
    public CMIComplete(int numPerm,ArrayList<List<Double>> data, boolean considerZeros, long seed) {
        this.numPerm=numPerm;
        this.origData=data;
        maxDeep=Integer.MAX_VALUE;
        this.considerZeros = considerZeros;
        if(seed != -1) random = new Random(seed);
        else random = ThreadLocalRandom.current();
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

        List<Double> toBeRandomized=origData.get(2);

        for (int i = 0; i < numPerm; i++) {
            List<Double> randomizedData= new ArrayList<>(toBeRandomized.size());
            ArrayList<Integer> currentPermutation = randomizedIndices.get(i);
            for (int j = 0; j < toBeRandomized.size(); j++) {
                randomizedData.add(toBeRandomized.get(currentPermutation.get(j)));
            }
            ArrayList<List<Double>> forPartitioning=new ArrayList<>(3);
            forPartitioning.add(origData.get(0));
            forPartitioning.add(origData.get(1));
            forPartitioning.add(randomizedData);
            IterativePartitioning ipRand=new IterativePartitioning(forPartitioning);
            randCMI[i]=ipRand.naivePartitioning();
        }
        int sum=0;
        for (int i = 0; i < numPerm; i++) {
            if(cmi<= randCMI[i]) sum++;
        }
        sum=Math.max(sum,1);
        pValue= (double) sum/numPerm;

    }

    public static double[] computeChiSquareCutoffs(double pSignificance, int dimensions){
        double[] chiSquareCutoffs = new double[dimensions];
        for(int dim = 0; dim < dimensions; dim++){
            double df = Math.pow(2, dim + 1) - 1;
            ChiSquaredDistribution chisquare = new ChiSquaredDistribution(df);
            double chi = chisquare.inverseCumulativeProbability(1 - pSignificance);
            chiSquareCutoffs[dim] = chi;
        }
        return(chiSquareCutoffs);
    }


    /**
     * This function should be used for CMI computation by IterativePartitioning and p-value estimation by randomization
     * of the third data array. This implementation uses the most efficient version of IterativePartitioning.
     */
    public void computeIterativePartitioning(Cube initialCube, double pSignificance){
        IterativePartitioning ip=new IterativePartitioning(origData);
        ip.setMaxDeep(maxDeep);
        ip.setConsiderZeros(considerZeros);
        ip.setChiSquareCutoffs(computeChiSquareCutoffs(pSignificance, origData.size()));
        cmi=ip.iterativePartitioningBetter(initialCube);

        long sum = IntStream.range(1, numPerm)
                .mapToDouble(i -> ComputeCMI.computeRandomCMI(origData.get(2).size(), ip, initialCube, random))
                .filter(randCMI -> randCMI >= cmi)
                .count();
        sum = Math.max(sum, 1);
        pValue = (double) sum/numPerm;
    }

    /**
     * This function uses CMI computation on uniform grid
     * @param numberOfBins number of bins of uniform grid in every dimension
     */
    public void computeUniformGrid(int numberOfBins){
        CMIUniform cmiUniform=new CMIUniform(origData);
        cmi=cmiUniform.computeCMI(numberOfBins);

        long sum = IntStream.range(1, numPerm).mapToDouble(
                i -> ComputeCMI.computeRandomCMIinUniformGrid(origData.get(2).size(), numberOfBins, origData, random))
                .filter(randCMI -> randCMI >= cmi)
                .count();
        sum = Math.max(sum, 1);
        pValue = (double) sum/numPerm;
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

        long sum = IntStream.range(1, numPerm).mapToDouble(
                i -> ComputeCMI.computeRandomCMIinPseudoUniformGrid(origData.get(2).size(), numberOfBins, ip, random))
                .filter(randCMI -> randCMI >= cmi)
                .count();
        sum = Math.max(sum, 1);
        pValue = (double) sum/numPerm;
        sum=Math.max(sum, 1);
        pValue=(double) sum/numPerm;
    }
}
