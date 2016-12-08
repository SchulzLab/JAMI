import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by fuksova on 12/1/15.
 * In this sparse matrix, score of interaction between genes and miRNA can be stored. Java sparse arrays are used
 * for storing the data. That is, rows do not have to have the same length!
 */
public class SparseScoreMatrix {

    double [][]values;  //Sparse array where double values are stored. Important: rows do not all have the same length
    int [][] indices; //Sparse array for storing column indices for every row. Access to values: M(i,j)=values[i][indices[i][j]]
    int [][] columnStructure; //Transposed information to indices. First index denotes column, second denotes rows for which the column contains nonzero values
    int maxRow; //number of rows in the complete matrix
    int maxColumn; //number of columns in the complete matrix


    /**
     * Initializes sparse matrix
     * @param indicesList list of matrix rows, stores indices with nonzero values.
     * @param valuesList list of matrix rows, stores values
     */
    public SparseScoreMatrix(ArrayList<int[]> indicesList,ArrayList<double[]> valuesList){ //For memory saving purposes, change to ArrayList of int[]/double[], I am affraid that this way memory could be duplicated
        values=new double[indicesList.size()][];                                            //Maybe it is worth to find out
        indices=new int[indicesList.size()][];

        maxColumn=-1;
        maxRow=valuesList.size();
        ArrayList<Integer> counter=new ArrayList<>();
        for (int i = 0; i < maxRow; i++) {
            values[i]=valuesList.get(i);
            int[] currentIndices = indicesList.get(i);
//            indices[i]=new int[ currentIndices.length];  //Does not work, it is necessary to use object Row
//            for (int j = 0; j < currentIndices.length; j++) {
//                indices[i][j]=currentIndices[j];
//            }
            indices[i]=currentIndices;

            int maxInd = currentIndices[currentIndices.length - 1];
            if(maxInd >maxColumn){
                for (int j = maxColumn+1; j <= maxInd; j++) {
                    counter.add(j,0);
                }
                maxColumn=maxInd;

            }
            for (int j = 0; j < currentIndices.length; j++) {
                int col=currentIndices[j];
                int act=counter.get(col);
                counter.set(col,++act);
            }
        }
        createColumnStructure(counter);
    }


    /**
     *
     * @param hmValues Hash map with score values stored according to gene and miRNA name created in ScoreStructure
     * @param intToGeneNames Map that express which row index in the matrix corresponds to certain gene name
     * @param intTomiRNANames Map that express which column index in the matrix corresponds to certain gene name
     */
    public SparseScoreMatrix(HashMap<String,HashMap<String,Double>> hmValues,ArrayList<String> intToGeneNames,ArrayList<String> intTomiRNANames){
        maxRow=intToGeneNames.size();
        maxColumn=intTomiRNANames.size();
        indices=new int[maxRow][];
        values=new double[maxRow][];
        ArrayList<Integer> counter=new ArrayList<>(maxColumn);
        for (int i = 0; i < maxColumn; i++) {
            counter.add(0);
        }
        for (int i = 0; i < maxRow; i++) {
            HashMap<String, Double> hmMiRNA = hmValues.get(intToGeneNames.get(i));
            indices[i]=new int[hmMiRNA.size()];
            values[i]=new double[hmMiRNA.size()];
            int columnCounter=0;
            for (int j = 0; j < maxColumn; j++) {
                Double value = hmMiRNA.get(intTomiRNANames.get(j));
                if(value!=null&&!value.equals(0.0)){
                    indices[i][columnCounter]=j;
                    values[i][columnCounter]=value;
                    columnCounter++;
                    int act=counter.get(j);
                    counter.set(j,++act);
                }
            }
        }
        createColumnStructure(counter);

    }


    /**
     * From the initialized sparse matrix structure creates column structure. That is, a list of nonzero rows for every
     * column
     * @param counter number of nonzero rows for given column
     */
    private void createColumnStructure(ArrayList<Integer> counter){
        columnStructure =new int[maxColumn+1][];
        for (int i = 0; i < maxColumn; i++) {
            columnStructure[i]=new int[counter.get(i)];
        }
        int [] counter2=new int[maxColumn+1];
        for (int i = 0; i < maxRow; i++) {
            int[] currentRow = indices[i];
            for (int j = 0; j < currentRow.length; j++) {
                int columnIndex = currentRow[j];
                columnStructure[columnIndex][counter2[columnIndex]]=i;
                counter2[columnIndex]++;
            }
        }
    }






    /**
     * Used for searching whether there is a nonzero value at given position.
     * @param rowIndex row index in full matrix
     * @param columnIndex column index in full matrix
     * @return returns column index in values where the value is stored or -1 if there is zero value
     */
    private int searchColumnValue(int rowIndex,int columnIndex){ //change to private
        int[] row=indices[rowIndex];
        int currentMaxColumn=row[row.length-1];
        if(columnIndex>currentMaxColumn){
            return -1;
        }
        int max=Math.min(columnIndex,row.length-1);
        int min=Math.max(0,row.length-1-currentMaxColumn+columnIndex);
        if(row[max]==columnIndex){
            return max;
        }
        if(row[min]==columnIndex){
            return min;
        }
        while(min<max){
            int center=(min+max)/2;
            if(row[center]==columnIndex){
                return center;
            }
            else if(row[center]>columnIndex){
                min=Math.max(min,center-row[center]+columnIndex);
                max=center-1;
            }
            else{
                max=Math.min(max,columnIndex-row[center]+center);
                min=center+1;
            }
        }
        if(min==max&&row[min]==columnIndex){
            return min;
        }
        else{
            return -1;
        }
    }

    /**
     *
     * @param rowIndex row index in full matrix
     * @param columnIndex column index in full matrix
     * @return value in the full matrix at given indices
     */
    public double getElementAt(int rowIndex,int columnIndex){
        int sparseColumnIndex=searchColumnValue(rowIndex,columnIndex);
        if(sparseColumnIndex>=0){
            return values[rowIndex][sparseColumnIndex];
        }
        else{
            return 0.0;
        }
    }

    /**
     *
     * @param rowIndex
     * @return indices of nonzero columns in given row
     */
    public int[] getIndicesAtRow(int rowIndex){
        return indices[rowIndex];
    }

    /**
     *
     * @param rowIndex
     * @return nonzero values in given row
     */
    public double[] getValuesAtRow(int rowIndex){
        return values[rowIndex];
    }
}
