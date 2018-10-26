package org.mpii.jami.helpers;

import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.mpii.jami.model.Triplet;
import java.util.List;

/**
 * Created by mlist on 11/3/17.
 */
public class BenjaminiHochberg {

    public static List<Triplet> adjustPValues(List<Triplet> triplets){

        double[] adjP = adjustPValues(triplets
                .stream()
                .sequential()
                .mapToDouble(t -> t.getpValue())
                .toArray());

        for(int i = 0; i < adjP.length; i++){
            triplets.get(i).setpAdjust(adjP[i]);
        }

        return(triplets);
    }

    public static double[] adjustPValues(double[] pValues){

        NaturalRanking ranking = new NaturalRanking(TiesStrategy.SEQUENTIAL);

        double[] pValueRanks = ranking.rank(pValues);
        int[] pValueIndices = new int[pValues.length];

        for(int i = 0; i < pValues.length; i++){
            pValueIndices[pValues.length - (int) pValueRanks[i]] = i;
        }
        double[] pValuesAdjusted = new double[pValues.length];

        double cumulativeMin = 1.0;

        for(int i = 0; i < pValues.length; i++){
            int k = pValueIndices[i];

            pValuesAdjusted[k] = Math.min(1.0, Math.min(cumulativeMin, pValues[k] * (pValues.length / pValueRanks[k])));
            cumulativeMin = Math.min(cumulativeMin, pValuesAdjusted[k]);
        }

        return(pValuesAdjusted);
    }
}
