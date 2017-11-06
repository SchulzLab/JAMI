package org.mpii.jami.helpers;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.util.List;

/**
 * Created by mlist on 11/6/17.
 */
public class FisherMethod {

    public static double combinedPValue(List<Double> pValues){
        ChiSquaredDistribution chi = new ChiSquaredDistribution(2 * pValues.size());
        double metaP = 1.0d - chi.cumulativeProbability(-2 * pValues.stream()
                .mapToDouble(i -> Math.log(i))
                .sum());

        return(metaP);
    }
}
