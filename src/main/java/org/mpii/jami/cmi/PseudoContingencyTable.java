package org.mpii.jami.cmi;

import org.apache.commons.math3.special.Beta;

/**
 * Created by fuksova on 11/23/15.
 * Currently only static methods are used, other functioins were replaced in ScoreStructure
 */
public class PseudoContingencyTable {
    double[][] table;

    public PseudoContingencyTable(double[][] scores){
        table=new double[2][2];
        double firstRow;
        double firstColumn;
        double sumAll;
        double firstElement;
        double[][] opposite = oppositeTable(scores);
        table[0][0]=Math.ceil(innerProduct(scores,scores));
        table[1][1]=Math.ceil(innerProduct(opposite,opposite));
        table[0][1]=Math.ceil(innerProduct(scores,opposite));
        table[1][0]=Math.ceil(innerProduct(opposite,scores));

    }

    private double[][] oppositeTable(double [][]scores){
        double [][]opposite=new double[scores.length][scores[0].length];
        for (int i = 0; i < scores.length; i++) {
            for (int j = 0; j < scores[0].length; j++) {
                opposite[i][j]=1-scores[i][j];
            }
        }
        return opposite;
    }

    private double innerProduct(double[][] first,double[][] second){
        double sum=0;
        for (int i = 0; i < first.length; i++) {
            sum+=first[i][0]*second[i][1];
        }
        return sum;
    }

    /**
     * Computes Fischer exact test p-value for given contingency table
     * @param table contingency table
     * @return p-value
     */
    public static double doTest(double [][]table){
        double firstRow=table[0][0]+table[0][1];
        double firstColumn=table[0][0]+table[1][0];
        double allPoints=firstRow+table[1][0]+table[1][1];
        return  test(table[0][0],firstColumn,firstRow,allPoints);
    }

    /**
     * Fischer exact test, instead of table use sums obtained from table
     * @param firstCell value in first table cell
     * @param firstColumn number of points in the first column
     * @param firstRow  number of points in the first row
     * @param allPoints number of all points in the table
     * @return p-value
     */
    public static double test(double firstCell, double firstColumn, double firstRow, double allPoints){
        double min=Math.max(0,firstColumn+firstRow-allPoints);
        double max=Math.min(firstRow,firstColumn);
        double sum=0.0;
        for (double i = firstCell; i <=max ; i++) {
            double log=lnBinomial(firstRow,i)+lnBinomial(allPoints-firstRow,firstColumn-i)-lnBinomial(allPoints,firstColumn);
            sum+=Math.exp(log);
        }
        return sum;
    }

    /**
     *
     * @param n
     * @param k
     * @return logarithm of binomial number
     */
    public static double lnBinomial(double n, double k){
        double q=n-k;
        if(q<=0){
            return 0;
        }
        q=-Math.log(q)-Beta.logBeta(k+1,q);
        return q;
    }


}
