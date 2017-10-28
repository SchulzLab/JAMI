package org.mpii.jami.tasks;

import org.mpii.jami.cmi.CMIUniform;

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
public class CMIUniformGridRecursiveTask extends RecursiveTask<Integer> {
    private int permutations;
    private double cmi;
    private int numOfSamples;
    private ArrayList<double[]> origData;
    private static final int THRESHOLD = 1000;
    private int numberOfBins;

    public CMIUniformGridRecursiveTask(int permutations,
                                       double cmi,
                                       ArrayList<double[]> origData,
                                       int numOfSamples,
                                       int numberOfBins
                            ) {
        this.permutations = permutations;
        this.cmi = cmi;
        this.numOfSamples = numOfSamples;
        this.origData = origData;
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

    private Collection<CMIUniformGridRecursiveTask> createSubtasks() {
        List<CMIUniformGridRecursiveTask> dividedTasks = new ArrayList<>();
        dividedTasks.add(new CMIUniformGridRecursiveTask(
                permutations / 2, cmi, origData, numOfSamples, numberOfBins));
        dividedTasks.add(new CMIUniformGridRecursiveTask(
                permutations / 2 + permutations % 2,
                cmi, origData, numOfSamples, numberOfBins));
        return dividedTasks;
    }

    private long processing(int permutations) {
        return(IntStream.range(1, permutations).mapToDouble(i -> computeRandomCMIinUniformGrid(numOfSamples, numberOfBins, origData))
                .filter(randCMI -> randCMI >= cmi)
                .count());
    }

     public static double computeRandomCMIinUniformGrid(int numOfSamples, int numberOfBins, ArrayList<double[]> origData){

        List<Integer> currentPermutation = IntStream.range(0, numOfSamples)
                .boxed().collect(Collectors.toList());
        Collections.shuffle(currentPermutation);

        double[] randomized= new double[numOfSamples];

        for (int j = 0; j < numOfSamples; j++) {
            randomized[j]= origData.get(2)[currentPermutation.get(j)];
        }
        ArrayList<double[]> randData=new ArrayList<>();
        randData.add(origData.get(0));
        randData.add(origData.get(1));
        randData.add(randomized);
        CMIUniform cmiRand=new CMIUniform(randData);
        return(cmiRand.computeCMI(numberOfBins));
    }

}
