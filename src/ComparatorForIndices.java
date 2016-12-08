import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by fuksova on 11/10/15.
 * This class is used in IterativePartitioning and in CompareDoubleArrays. It is used for replacing the original
 * double data with their integer ascending orders separately in every dimension.
 */
public class ComparatorForIndices implements Comparator<Integer>{
    private double [] sourceArray;

    public ComparatorForIndices(double [] sourceArray){
        this.sourceArray=sourceArray;

    }


    public ComparatorForIndices(Double [] sourceArray){
        this.sourceArray=new double[sourceArray.length];
        for (int i = 0; i < sourceArray.length; i++) {
            this.sourceArray[i]=sourceArray[i];
        }
    }

    /**
     * Used to get array to be sorted using this comparator.
     * @return array of indices of the source array, that is from 0 to sourceArray.length-1
     */
    public Integer[] getArrayToSort(){
        Integer[] indices=new Integer[sourceArray.length];
        for (int i = 0; i < sourceArray.length; i++) {
            indices[i]=i;
        }
        return indices;
    }

    /**
     *
     * @param sortedIndices array of sorted indexes
     * @return inverse array to the input array
     */
    public Integer[] getInverseArray(Integer[] sortedIndices){
        Integer [] inverseArray=new Integer[sourceArray.length];
        for (int i = 0; i < sourceArray.length; i++) {
            inverseArray[sortedIndices[i]]=i;
        }
        return inverseArray;
    }

    /**
     * This way, the array with integer orders replacing the original double values is obtained directly.
     * @return array with integer values
     */
    public Integer[] computeInverseSorted(){
        Integer[] array = getArrayToSort();
        Arrays.sort(array,this);
        return getInverseArray(array);
    }

    @Override
    /**
     * This comparator is to be used exclusively for sorting indexes of the array used in constructor. The best way
     * to use it is sorting the array returned by getArrayToSort.
     */
    public int compare(Integer i1, Integer i2) {
        if(i1<0||i1>=sourceArray.length||i2<0||i2>=sourceArray.length){
            throw new IllegalArgumentException("Index out of bound. This comparator is to be used exclusively for sorting indexes of the array used in constructor.");
        }
        else{
            return Double.compare(sourceArray[i1],sourceArray[i2]);
        }
    }
}
