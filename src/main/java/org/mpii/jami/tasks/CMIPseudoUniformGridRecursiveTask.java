package org.mpii.jami.tasks;

import org.mpii.jami.cmi.IterativePartitioning;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveTask;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by mlist on 12/13/16.
 */
public class CMIPseudoUniformGridRecursiveTask extends RecursiveTask<Integer> {
    private int permutations;
    private double cmi;
    private int numOfSamples;
    private IterativePartitioning ip;
    private static final int THRESHOLD = 1000;
    private int numberOfBins;

    public CMIPseudoUniformGridRecursiveTask(int permutations,
                                             double cmi,
                                             IterativePartitioning ip,
                                             int numOfSamples,
                                             int numberOfBins
                            ) {
        this.permutations = permutations;
        this.cmi = cmi;
        this.numOfSamples = numOfSamples;
        this.ip = ip;
        this.numberOfBins = numberOfBins;
    }

    @Override
    protected Integer compute() {
        if (permutations > THRESHOLD) {
            return ForkJoinTask.invokeAll(createSubtasks())
            .stream().mapToInt(ForkJoinTask::join)
            .sum();
        } else {
            return((int) processing(permutations));
        }
    }

    private Collection<CMIPseudoUniformGridRecursiveTask> createSubtasks() {
        List<CMIPseudoUniformGridRecursiveTask> dividedTasks = new ArrayList<>();
        dividedTasks.add(new CMIPseudoUniformGridRecursiveTask(
                permutations / 2, cmi, ip, numOfSamples, numberOfBins));
        dividedTasks.add(new CMIPseudoUniformGridRecursiveTask(
                permutations / 2 + permutations % 2,
                cmi, ip, numOfSamples, numberOfBins));
        return dividedTasks;
    }

    private long processing(int permutations) {
        return(IntStream.range(1, permutations).mapToDouble(i -> computeRandomCMIinPseudoUniformGrid(numOfSamples, numberOfBins, ip))
                .filter(randCMI -> randCMI >= cmi)
                .count());
    }

     public static double computeRandomCMIinPseudoUniformGrid(int numOfSamples, int numberOfBins, IterativePartitioning ip){

        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation);

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

}
