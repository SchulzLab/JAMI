import java.util.Random;

/**
 * Created by fuksova on 1/20/16.
 * This class can be used for comparison of two double arrays to find out whether the orderings of the two arrays is
 * similar or not. For that purpose several methods are provided.
 */
public class CompareDoubleArrays {

    double kendall; //Value of Kendall rank correlation coefficient
    double avgDifference;  //Average absolute value of difference between entries of the two vectors
    double avgOrderDifference;  //Average absolute value difference between integer orders of the two vectors
    double minDifference;
    int minOrderDifference;
    double maxDifference;
    int maxOrderDifference;
    double [] data1;
    double [] data2;
    Integer[] order1;
    Integer[] order2;

    public CompareDoubleArrays(double[] data1, double[] data2){
        if(data1.length!=data2.length){
            throw new IllegalArgumentException("data1 and data2 must have the same length");
        }
        this.data1=data1;
        this.data2=data2;
    }

    public static double kendall(double[] firstArray,double[] secondArray){
        int n=firstArray.length;
        double stat=0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j <= i - 1; j++) {
                stat=stat+Math.signum((firstArray[i]-firstArray[j]))*Math.signum((secondArray[i] - secondArray[j]));
            }
        }
        return stat*2.0/(n*(n-1));
    }

    public void computeAll(){
        avgDiff();
        avgOrderDiff();
        minDifference();
        maxDifference();
        minOrderDifference();
        maxOrderDifference();
        kendall=kendall(data1,data2);
    }

    public double avgOrderDiff(){
        initOrders();
        double sum=0.0;
        for (int i = 0; i < order1.length; i++) {
            sum+=Math.abs(order1[i]-order2[i]);
        }
        avgOrderDifference=sum/(double)order1.length;
        return avgDifference;
    }

    public double avgDiff(){
        double sum=0.0;
        for (int i = 0; i < data1.length; i++) {
           sum+=Math.abs(data1[i]-data2[i]);
        }
        avgDifference=sum/(double)data1.length;
        return avgDifference;
    }

    public int minOrderDifference(){
        initOrders();
        int min=order1.length;
        for (int i = 0; i < order1.length; i++) {
            int temp=Math.abs(order1[i]-order2[i]);
            if(temp<min){
                min=temp;
            }
        }
        minOrderDifference=min;
        return minOrderDifference;
    }

    public int maxOrderDifference(){
        initOrders();
        int max=0;
        for (int i = 0; i < order1.length; i++) {
            int temp=Math.abs(order1[i]-order2[i]);
            if(temp>max){
                max=temp;
            }
        }
        maxOrderDifference=max;
        return maxOrderDifference;
    }

    public double maxDifference(){
        double max=Double.MIN_VALUE;
        for (int i = 0; i < data1.length; i++) {
            double temp=Math.abs(data1[i]-data2[i]);
            if(temp>max){
                max=temp;
            }
        }
        maxDifference=max;
        return maxDifference;
    }

    public double minDifference(){
        double min =Double.MAX_VALUE;
        for (int i = 0; i < data1.length; i++) {
            double temp=Math.abs(data1[i]-data2[i]);
            if(temp<min){
                min=temp;
            }
        }
        minDifference=min;
        return minDifference;
    }

    /**
     * Computes all values and outputs them
     */
    public void computeAndOutput(){
        computeAll();

        System.out.println("avg difference: "+avgDifference);
        System.out.println("avg order difference: "+avgOrderDifference);
        System.out.println("max difference: "+maxDifference);
        System.out.println("min difference: "+minDifference);
        System.out.println("max order difference: "+maxOrderDifference);
        System.out.println("min order difference: "+minOrderDifference);
        System.out.println("kendall: "+kendall);
    }


    /**
     * Transform double values to integer orders
     */
    public void initOrders(){
        if(order1==null){
            ComparatorForIndices ci=new ComparatorForIndices(data1);
            order1=ci.computeInverseSorted();
        }
        if(order2==null){
            ComparatorForIndices ci=new ComparatorForIndices(data2);
            order2=ci.computeInverseSorted();
        }
    }

    public static double[] randArray(int lengthOfArray){
        double [] randA=new double[lengthOfArray];
        Random random=new Random();
        for (int i = 0; i < lengthOfArray; i++) {
            //randA[i]=random.nextDouble();
            randA[i]=Math.max(random.nextInt(1001),1.0)/1000.0;
        }
        return randA;
    }

     /**
     * How the results look like for random arrays
     * @param numberOfRuns how many random arrays
     * @param lengthOfArray how long arrays
     */
    public static void tryTwoArrayStatistics(int numberOfRuns,int lengthOfArray){
        double[] kendall=new double[numberOfRuns];
        double[] avgDifference=new double[numberOfRuns];
        double[] avgOrderDifference=new double[numberOfRuns];
        for (int i = 0; i < numberOfRuns; i++) {
            double []array1=randArray(lengthOfArray);
            double []array2=randArray(lengthOfArray);
            CompareDoubleArrays evaluate=new CompareDoubleArrays(array1,array2);
            evaluate.computeAll();
            kendall[i]=evaluate.kendall;
            avgDifference[i]=evaluate.avgDifference;
            avgOrderDifference[i]=evaluate.avgOrderDifference;
        }
        double meanK=AdditionalComputations.computeMean(kendall);
        double stdKendal=AdditionalComputations.computeStdDev(meanK,kendall);
        double meanAvgDiff=AdditionalComputations.computeMean(avgDifference);
        double stdAvgDiff=AdditionalComputations.computeStdDev(meanAvgDiff,avgDifference);
        double meanAvgOrderDiff=AdditionalComputations.computeMean(avgOrderDifference);
        double stdAvgOrderDiff=AdditionalComputations.computeStdDev(meanAvgOrderDiff,avgOrderDifference);
        System.out.println("Kendall: mean, st. dev: "+meanK+", "+stdKendal);
        System.out.println("Avg. diff.: mean, st. dev: "+meanAvgDiff+", "+stdAvgDiff);
        System.out.println("Avg. order diff.: mean, st. dev: "+meanAvgOrderDiff+", "+stdAvgOrderDiff);


    }


}
