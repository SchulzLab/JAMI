package org.mpii.jami.helpers;

/**
 * Created by fuksova on 1/22/16.
 * Some basic mathematical methods that can be useful
 */
public class AdditionalComputations {

    /**
     * Computes mean
     * @param data input data
     * @return
     */
    public static double computeMean(double[] data){
        double sum=0.0;
        for (int i = 0; i < data.length; i++) {
            sum+=data[i];
        }
        return sum/(double)data.length;
    }

    /**
     * Standard deviation with precomputed mean
     * @param mean
     * @param data
     * @return
     */
    public static double computeStdDev(double mean,double[] data){
        double sum=0.0;
        for (int i = 0; i < data.length; i++) {
            sum+=(data[i]-mean)*(data[i]-mean);
        }
        sum=sum/data.length;
        return Math.sqrt(sum);
    }

    /**
     * Standard deviation without precomputed mean
     * @param data
     * @return
     */
    public static double computeStdDev(double[] data){
        double mean=computeMean(data);
        return computeStdDev(mean,data);
    }
}
