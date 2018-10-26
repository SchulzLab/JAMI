package org.mpii.jami.helpers;

import org.mpii.jami.cmi.CMIUniform;
import org.mpii.jami.cmi.Cube;
import org.mpii.jami.cmi.IterativePartitioning;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ComputeCMI {

    public static double computeRandomCMI(int numOfSamples, IterativePartitioning ip, Cube initialCube, Random random) {
        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation, random);

        Integer[] randomizedSorted = new Integer[numOfSamples];
        Integer[] randomizedInverse = new Integer[numOfSamples];

        for (int j = 0; j < numOfSamples; j++) {
            Integer value = ip.getInverseSortedIndices().get(2)[currentPermutation.get(j)];
            randomizedInverse[j] = value;
            randomizedSorted[value] = j;
        }

        ArrayList<Integer[]> sorted = new ArrayList<>();
        sorted.add(ip.getSortedIndices().get(0));
        sorted.add(ip.getSortedIndices().get(1));
        sorted.add(randomizedSorted);
        ArrayList<Integer[]> inverse = new ArrayList<>();
        inverse.add(ip.getInverseSortedIndices().get(0));
        inverse.add(ip.getInverseSortedIndices().get(1));
        inverse.add(randomizedInverse);
        IterativePartitioning ipRand = new IterativePartitioning(sorted, inverse);
        ipRand.setMaxDeep(ip.getMaxDeep());
        ipRand.setConsiderZeros(ip.ConsiderZeros());
        ipRand.setMinimumValues(ip.getMinimumValues());
        ipRand.setChiSquareCutoffs(ip.getChiSquareCutoffs());
        return (ipRand.iterativePartitioningBetter(initialCube));
    }

    public static double computeRandomCMIinPseudoUniformGrid(int numOfSamples, int numberOfBins, IterativePartitioning ip, Random random){

        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation, random);

        Integer[] randomizedSorted=new Integer[numOfSamples];
        Integer[] randomizedInverse=new Integer[numOfSamples];
        for (int j = 0; j < numOfSamples; j++) {
            Integer value = ip.getInverseSortedIndices().get(2)[currentPermutation.get(j)];
            randomizedInverse[j]= value;
            randomizedSorted[value]=j;
        }

        ArrayList<Integer[]> sorted=new ArrayList<>();
        sorted.add(ip.getSortedIndices().get(0));
        sorted.add(ip.getSortedIndices().get(1));
        sorted.add(randomizedSorted);
        ArrayList<Integer[]> inverse=new ArrayList<>();
        inverse.add(ip.getInverseSortedIndices().get(0));
        inverse.add(ip.getInverseSortedIndices().get(1));
        inverse.add(randomizedInverse);
        IterativePartitioning ipRand=new IterativePartitioning(sorted,inverse);

        return(ipRand.cmiInUniformGrid(numberOfBins));
    }

    public static double computeRandomCMIinUniformGrid(int numOfSamples, int numberOfBins, ArrayList<List<Double>> origData, Random random){

        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation, random);

        List<Double> randomized= new ArrayList<>(numOfSamples);

        for (int j = 0; j < numOfSamples; j++) {
            randomized.add(origData.get(2).get(currentPermutation.get(j)));
        }
        ArrayList<List<Double>> randData=new ArrayList<>();
        randData.add(origData.get(0));
        randData.add(origData.get(1));
        randData.add(randomized);
        CMIUniform cmiRand=new CMIUniform(randData);
        return(cmiRand.computeCMI(numberOfBins));

    }
}
