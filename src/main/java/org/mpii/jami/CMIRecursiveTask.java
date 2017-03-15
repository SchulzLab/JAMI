package org.mpii.jami;

import java.util.*;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveTask;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by mlist on 12/13/16.
 */
public class CMIRecursiveTask extends RecursiveTask<Integer> {
    private int permutations;
    private double cmi;
    private int numOfSamples;
    private IterativePartitioning ip;
    private static final int THRESHOLD = 100;

    public CMIRecursiveTask(int permutations,
                            double cmi,
                            IterativePartitioning ip,
                            int numOfSamples
                            ) {
        this.permutations = permutations;
        this.cmi = cmi;
        this.numOfSamples = numOfSamples;
        this.ip = ip;
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

    private Collection<CMIRecursiveTask> createSubtasks() {
        List<CMIRecursiveTask> dividedTasks = new ArrayList<>();
        dividedTasks.add(new CMIRecursiveTask(
                permutations / 2, cmi, ip, numOfSamples));
        dividedTasks.add(new CMIRecursiveTask(
                permutations / 2 + permutations % 2,
                cmi, ip, numOfSamples));
        return dividedTasks;
    }

    private long processing(int permutations) {
        return(IntStream.range(1, permutations).mapToDouble(i -> computeRandomCMI())
                .filter(randCMI -> randCMI <= cmi)
                .count());
    }

    private double computeRandomCMI() {
        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation);

        Integer[] randomizedSorted = new Integer[numOfSamples];
        Integer[] randomizedInverse = new Integer[numOfSamples];

        for (int j = 0; j < numOfSamples; j++) {
            Integer value = ip.inverseSortedIndices.get(2)[currentPermutation.get(j)];
            randomizedInverse[j] = value;
            randomizedSorted[value] = j;
        }

        ArrayList<Integer[]> sorted = new ArrayList<>();
        sorted.add(ip.sortedIndices.get(0));
        sorted.add(ip.sortedIndices.get(1));
        sorted.add(randomizedSorted);
        ArrayList<Integer[]> inverese = new ArrayList<>();
        inverese.add(ip.inverseSortedIndices.get(0));
        inverese.add(ip.inverseSortedIndices.get(1));
        inverese.add(randomizedInverse);
        IterativePartitioning ipRand = new IterativePartitioning(sorted, inverese);
        ipRand.maxDeep = ip.maxDeep;
        return (ipRand.iterativePartitioningBetter());
    }

}