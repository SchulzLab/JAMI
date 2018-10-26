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

        NaturalRanking ranking = new NaturalRanking(TiesStrategy.MAXIMUM);

        double[] pValueRanks = ranking.rank(pValues);
        double[] pValuesAdjusted = new double[pValues.length];

        for(int i = 0; i < pValues.length; i++){
            pValuesAdjusted[i] = Math.min(pValues[i] * (pValues.length / pValueRanks[i]),
                    pValues[i+1] * (pValues.length / pValueRanks[i+1]));
        }

        return(pValuesAdjusted);
    }
}
