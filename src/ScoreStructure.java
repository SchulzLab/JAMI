import java.io.*;
import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by fuksova on 12/9/15.
 * This class contains methods for reading Cupid score files, computing p-values and q-values for all gene pairs and
 * for obtaining p-value threshold in the same way that is used in Cupid. However, all those methods should be used
 * carefully because there is no warranty about the score correctness.
 * Other problems: The methods used for p-value computation are quite strange heuristics. They prefer those gene pairs
 * that have very low score value with almost all miRNA. Moreover, their threshold chooses over 3 millions gene pairs
 * to be further tested.
 *
 * Q-values computation is based on software used in
 * Kall et al. 2008: Non-parametric estimation of posterior error probabilities associated with peptides identified by
 * tandem mass spectrometry
 * The used method is not described in the paper, I implemented it based on the source code of the software.
 */
public class ScoreStructure {

    /*


    Klidne by tu mohl byt i seznam vice souboru s daty, mozna by se prislo i na par dalsich statickych promennych,
    ktere by tu mohly byt... - treba separator (nebo i nestatickych, ted uz by to mohlo byt jedno...?)
    readFiles by nemuselo
     */


    String positiveFilename;  //File name where score for positive examples can be found
    String negativeFilename;  //File name where score for other examples can be found
    public static int geneInd;  //Index of column containing gene names
    public static int miRNAInd;  //Index of column containing miRNA names
    public static int valueInd;  //Index of column containing score values
    private ArrayList<String> intToGeneName; //Mapping of integers to corresponding gene names (for identifying a suitable row index in scoreMatrix)
    private ArrayList<String> intTomiRNAName; //Mapping of integers to corresponding miRNA names (for identifying a suitable column index in scoreMatrix)
    HashMap<String,Integer> geneNameToInt; //Opposite map to intToGeneName
    HashMap<String,Integer> miRNANameToInt; //Opposite map to intTomiRNAName
    SparseScoreMatrix scoreMatrix; //Matrix with Cupid step 2 score. Rows are genes, columns are miRNAs. The above maps are used for identifying indices for corresponding names
    int genesNumber; //Number of genes in score matrix
    public static Random rand=new Random(); //for randomized operations

    //public static String separator;


    /**
     *
     * @param filenamePositive name of file with positive examples
     * @param filenameNegative name of fiel with other examples
     */
    public ScoreStructure(String filenamePositive,String filenameNegative){

        this.negativeFilename=filenameNegative;
        this.positiveFilename=filenamePositive;

    }




    /**
     *Initialization of scoreMatrix and name to integer mapping fields
     * @param skipFirstRowPos If true first row in the file with positive examples will be skipped
     * @param skipFirstRowNeg If true first row in the file with other examples will be skipped
     */
    public void init(boolean skipFirstRowPos,boolean skipFirstRowNeg){
        long start=System.currentTimeMillis();
        HashMap<String, HashMap<String, Double>> hm = readPosAndNegData(positiveFilename, negativeFilename,skipFirstRowPos,skipFirstRowNeg);
        long end=System.currentTimeMillis();
        System.out.println("Time: "+(end-start));
        intToGeneName=new ArrayList<>();
        HashSet<String> miRNANamesSet=new HashSet<>();
        for (Map.Entry<String, HashMap<String, Double>> geneEntry : hm.entrySet()) {
            String geneName=geneEntry.getKey();
            HashMap<String, Double> miRNAMap = geneEntry.getValue();
            intToGeneName.add(geneName);
            for (String miRNAName : miRNAMap.keySet()) {
                miRNANamesSet.add(miRNAName);
            }
        }
        intTomiRNAName=new ArrayList<>(miRNANamesSet);
        Collections.sort(intTomiRNAName);
        Collections.sort(intToGeneName);
        initHashMaps();

        scoreMatrix=new SparseScoreMatrix(hm,intToGeneName,intTomiRNAName);


    }

    /**
     * Used in init()
     */
    private void initHashMaps(){
        geneNameToInt=new HashMap<>();
        miRNANameToInt=new HashMap<>();
        genesNumber=intToGeneName.size();
        for (int i = 0; i < genesNumber; i++) {
            geneNameToInt.put(intToGeneName.get(i),i);
        }
        for (int i = 0; i < intTomiRNAName.size(); i++) {
            miRNANameToInt.put(intTomiRNAName.get(i),i);
        }
    }




    /**
     * If there is initialized hash map with score values read from one file, this method enables to append score
     * data from other file to it.
     * @param allOnes Instead of reading exact score value from file, use 1.0
     * @param fileName Name of file
     * @param oldData HashMap initialized previously from other score file or is empty
     * @param skipFirstRow Indicates whether the first row of the file should be skipped
     */
    public static void appendScoreData(boolean allOnes, String fileName, HashMap<String, HashMap<String, Double>> oldData,boolean skipFirstRow) {
        Pattern pattern = Pattern.compile("\t");   //should be possible to set. Maybe a couple of static variables in this class? Instead of this big number of input parameters!
        FileReader fr;
        try {
            fr = new FileReader(fileName);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
            if(skipFirstRow){
                line=br.readLine();
            }
            while (line != null) {
                String[] entries = pattern.split(line);
                String geneID=entries[geneInd];
                String miRNAID=entries[miRNAInd];
//                if(miRNAID.equals("hsa-miR-449b-5p")){
//                    System.out.println(miRNAID);
//                }
                Double value=1.0;
                if(!allOnes){
                    value=Double.parseDouble(entries[valueInd]);
                }
                HashMap<String, Double> map = oldData.get(geneID);
                if(map==null){
                    map=new HashMap<>();
                    oldData.put(geneID,map);
                }
                else if(map.containsKey(miRNAID)){
                    double oldValue=map.get(miRNAID);
                    value=Math.max(oldValue,value);
                }
                if(!value.equals(0.0)) {
                    map.put(miRNAID, value);
                }
                line=br.readLine();
            }

            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Reads score data from a file and returns them in a form of HashMap
     * @param allOnes Instead of values from the score file 1.0 is used as default
     * @param fileName Name of file
     * @param skipFirstRow Skip first row
     * @return HashMap with score data geneID->miRNAID->value
     */
    public static HashMap<String,HashMap<String,Double>> readScoreData(boolean allOnes, String fileName,boolean skipFirstRow){
        HashMap<String,HashMap<String,Double>> hm=new HashMap<>();
        appendScoreData(allOnes, fileName, hm,skipFirstRow);
        return hm;
    }

    /**
     * Reads score data from two files, the one contains positive examples and the other contains other examples
     * @param posFileName file name with positive examples
     * @param negFileName file name with other examples
     * @param skipFirstRowPos Skip first row in file with positive examples
     * @param skipFirstRowNeg Skip first row in file with other examples
     * @return HashMap with score data in a form <gene id,<miRNAid,score value>>
     */
    public static HashMap<String,HashMap<String,Double>> readPosAndNegData(String posFileName,String negFileName,boolean skipFirstRowPos,boolean skipFirstRowNeg){ //combines the previous two functions
        HashMap<String, HashMap<String, Double>> hm = readScoreData(true, posFileName,skipFirstRowPos);
        appendScoreData(false,negFileName,hm,skipFirstRowNeg);
        return hm;

    }


    /**
     * Computes pseudo contingency table as is used in Cupid
     * @param gene1 index of first gene
     * @param gene2 index of second gene
     * @return pseudo contingency table
     */
    public double[][] computePseudoContTable(int gene1, int gene2){
        if(gene1==gene2||gene1>=scoreMatrix.maxRow||gene2>= scoreMatrix.maxRow){
            throw new IllegalArgumentException("Wrong combination of indices");
        }
        double[][] table=new double[2][2];
        int [] row1Indices=scoreMatrix.getIndicesAtRow(gene1);
        double [] row1Values=scoreMatrix.getValuesAtRow(gene1);
        int [] row2Indices=scoreMatrix.getIndicesAtRow(gene2);
        double [] row2Values=scoreMatrix.getValuesAtRow(gene2);
        int pointer1=0;
        int pointer2=0;
        for (int i = 0; i < intTomiRNAName.size(); i++) {
            if(pointer1>=row1Indices.length||i<row1Indices[pointer1]){
                if(pointer2>=row2Indices.length||i<row2Indices[pointer2]){ //will cause problems when pointer is out of bounds -> extract variable...
                    table[1][1]++;
                }
                else{
                    double value2=row2Values[pointer2];
                    pointer2++;
                    table[1][0]+=value2;
                    table[1][1]+=(1.0-value2);
                }
            }
            else{
                double value1=row1Values[pointer1];
                pointer1++;
                if(pointer2>=row2Indices.length||i<row2Indices[pointer2]){
                    table[1][1]+=(1-value1);
                    table[0][1]+=value1;

                }
                else{
                    double value2=row2Values[pointer2];
                    pointer2++;
                    table[0][0]+=value2*value1;
                    table[0][1]+=value1*(1.0-value2);
                    table[1][0]+=(1-value1)*value2;
                    table[1][1]+=(1-value1)*(1-value2);
                }
            }
        }
        table[0][0]=Math.ceil(table[0][0]);
        table[0][1]=Math.ceil(table[0][1]);
        table[1][0]=Math.ceil(table[1][0]);
        table[1][1]=Math.ceil(table[1][1]);
        return table;
    }

    /**
     * Computes p-value for the interaction between two genes based on pseudo contingency table. Please be aware that
     * this method is probably not suitable for gene pairs filtering.
     * @param gene1 first gene index
     * @param gene2 second gene index
     * @return p-value
     */
    public double computePValue(int gene1,int gene2){
        return PseudoContingencyTable.doTest(computePseudoContTable(gene1,gene2));
    }

    /**
     * Computes p-value for the interaction between two genes based on pseudo contingency table. Please be aware that
     * this way is probably not suitable for gene pairs filtering.
     * @param gene1Name name of first gene
     * @param gene2Name name of second gene
     * @return p-value
     */
    public double computePValue(String gene1Name,String gene2Name){
        return computePValue(geneNameToInt.get(gene1Name),geneNameToInt.get(gene2Name));
    }


    /**
     * Computes all p-values and return them. This procedure is really time consuming for large number of gene pairs
     * @return array of p-values in certain order, procedure getIndicesToOrder identifies the gene pair corresponding
     * to certain index in this returned array
     */
    public double [] getAllPValues(){
        int sizeOfArray=genesNumber*(genesNumber-1)/2;
        double [] pValues=new double[sizeOfArray];
        int counter=0;
        for (int i = 0; i < genesNumber -1; i++) {
            for (int j = i+1; j < genesNumber; j++) {
                pValues[counter]=computePValue(i,j);
                int[] check=getIndicesToOrder(counter);
                if(check[0]!=i||check[1]!=j){
                   System.out.println(counter+" should be: "+i+" "+j+", but is "+ check[0]+", "+check[1]);
                }
                counter++;
            }
        }
        return pValues;
    }

    /**
     * p-values for all gene pairs are computed an stored in a file
     * @param fileName name of the file
     */
    public void writePValuesToFile(String fileName){
        FileWriter fw = null;
        double[] pValues=getAllPValues();
        try {
            fw = new FileWriter(fileName);
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < pValues.length; i++) {
                bw.write(pValues[i]+"\n");
                bw.flush();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * P-values that are given as input fo this function are stored in a file
     * @param fileName name of the file
     * @param pValues computed p-values
     */
    public static void writePValuesToFile(String fileName, double[] pValues){
        FileWriter fw = null;
        try {
            fw = new FileWriter(fileName);
            BufferedWriter bw = new BufferedWriter(fw);
            for (int i = 0; i < pValues.length; i++) {
                bw.write(Double.toString(pValues[i]));
                bw.write("\n");
                bw.flush();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Based on p-values and threshold, this method selects gene pairs and stores them into a file
     * @param fileName name of the file where pairs should be stored
     * @param pValues computed p-values in a correct order regarding list of genes!
     * @param threshold p-values below this threshold are considered significant
     */
    public void pairsBellowThresholdToFile(String fileName, double[] pValues,double threshold){
        FileWriter fw = null;
        try {
            fw = new FileWriter(fileName);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("threshold: "+threshold+"\n");
            for (int i = 0; i < pValues.length; i++) {
                if(pValues[i]<=threshold){
                    int [] indices=getIndicesToOrder(i);
                    bw.write(intToGeneName.get(indices[0])+"\t"+intToGeneName.get(indices[1])+"\t"+pValues[i]);
                    bw.write("\n");
                    bw.flush();

                }
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Instead of computing p-values directly, reads them from a file. The file has to be in the format as the output file
     * of functions writePValuesToFile. If we want to use them together with correct gene names, it is necessary that gene names
     * are read from the same file as by creating the p-value file. This is the only way, how p-values can be correctly
     * reassigned.
     * @param fileName Name of file where p-values are stored
     * @param size number of lines (numerical values) in the file
     * @return
     */
    public static double[] getPValuesFromFile(String fileName,int size){
        double [] pValues=new double[size];
        FileReader fr;
        try {
            fr = new FileReader(fileName);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
            int counter=0;
            while(counter<size&&line!=null){
                //pValues[counter]=Double.parseDouble(line);  //Probably not right but I am also probably but I obtained pvalues larger than one
                pValues[counter]=Math.min(1.0, Double.parseDouble(line));  //Probably not right but I am also probably but I obtained pvalues larger than one
                counter++;
                line=br.readLine();
            }
            br.close();
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return pValues;
    }


    /**
     * In case that we have a set of p-values and a threshold on q-values, this method computes q-values from the list
     * of p-values and outputs corresponding p-value threshold.
     * @param qThreshold threshold on q-values
     * @param pValues list of p-values
     * @param filename If this is not null, sorted lists of p-values and corresponding q-values are stored into this file
     *                 They are stored in sorted orders. This file can not serve as source file for getPValuesFromFile!
     * @return p-value threshold
     */
    public static double getPThresholdToQThreshold(double qThreshold, double[] pValues, String filename) {
        double[] qvalues = getQvaluesFromPvalues(pValues, false, filename);
        int counter = 0;
        double pThreshold = 0.0;
        while (counter < qvalues.length) {
            if (qvalues[counter] >= qThreshold) {
                break;
            }
            counter++;
        }

        if (counter != 0) {
            pThreshold = pValues[counter - 1];
        }
        return pThreshold;
    }

    /**
     * Computes q-values based on a list of p-values
     * @param pValues list of p-values
     * @param preserveOriginal if true, the original list of p-values is preserved, otherwise is replaced with its sorted
     *                         version
     *                         If it is replaced with sorted version, original gene pairs can not be recovered for
     *                         certain p-value index.
     * @return list of q-values
     */
    public static double[] getQvaluesFromPvalues(double [] pValues,boolean preserveOriginal){
        return  getQvaluesFromPvalues(pValues,preserveOriginal,null);
    }


    /**
     * Computes q-values based on a list of p-values
     * @param pValues list of p-values
     * @param preserveOriginal if true, the original list of p-values is preserved, otherwise is replaced with its sorted version
     *                         If it is replaced with sorted version, original gene pairs can not be recovered for
     *                         certain p-value index.
     * @param filename If this is not null, sorted lists of p-values and corresponding q-values are stored into this file
     *                 They are stored in sorted orders. This file can not serve as source file for getPValuesFromFile!
     * @return list of q-values
     */
    public static double[] getQvaluesFromPvalues(double [] pValues,boolean preserveOriginal,String filename){
        double [] pValuesSorted;
        if(preserveOriginal){
            pValuesSorted=Arrays.copyOf(pValues,pValues.length);
           // Arrays.sort(pValuesSorted);
        }
        else{
            pValuesSorted=pValues;
        }
        Arrays.sort(pValuesSorted);

        int N = pValues.length;
        double [] qvalues=new double[N];
        double pi0=getPi0(pValuesSorted);
        System.out.println("pi0"+pi0);
        for (int i = 1; i <= N; i++) {
            qvalues[i-1]=pValuesSorted[i-1]* N *pi0/((double)(i));
        }
        double value=qvalues[N-1];
        for (int i = N-2; i >= 0; i--) {
            if(value<qvalues[i]){
                qvalues[i]=value;
            }
            else{
                value=qvalues[i];
            }
        }
        if(filename!=null){
            FileWriter fw = null;

            try {
                fw = new FileWriter(filename);
                BufferedWriter bw = new BufferedWriter(fw);
                bw.write("p-values\tq-values\n");
                for (int i = 0; i < pValues.length; i++) {
                    bw.write(Double.toString(pValuesSorted[i]));
                    bw.write("\t");
                    bw.write(Double.toString(qvalues[i]));
                    bw.write("\n");
                    bw.flush();
                }
                bw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return qvalues;
    }


    /**
     * This value is used for q-value computation. I do not know what is it supposed to be, this is just reimplemented
     * from original code from Kall et al.
     * @param pvalues list of p-values
     * @return pi0 value
     */
    private static double getPi0(double[]pvalues){
        int numBoot=100;
        double maxLambda=0.5;
        int numLambda=100;
        ArrayList<Double> pi0s,lambdas,mse,pBoot;
        pi0s=new ArrayList<>();
        lambdas=new ArrayList<>();
        mse=new ArrayList<>();
       // pBoot=new ArrayList<>();

        int start=0;
        for (int ix = 0; ix <= numLambda; ix++) {
            double lambda=((ix+1)/(double)numLambda)*maxLambda;
            start=lowerBound(pvalues,0,lambda);
            double wl=(double)(pvalues.length-start);
            double pi0=wl/pvalues.length/(1-lambda);
            if(pi0>0.0){
                lambdas.add(lambda);
                pi0s.add(pi0);
            }
        }
        if(pi0s.size()==0){
            System.out.println("Error in the input data: too good separation between targetand decoy. Impossible to estimate pi0. Terminating.");
            return -1;
        }
        double minPi0= Collections.min(pi0s);
        for (int i = 0; i < lambdas.size(); i++) {
            mse.add(0.0);
        }
        for (int j = 0; j < numBoot; j++) {
            pBoot=bootstrap(pvalues,1000);
            int n=pBoot.size();
            for (int i = 0; i < lambdas.size(); i++) {
                start=lowerBound(pBoot,0,lambdas.get(i));
                double wl=pBoot.size()-start;
                double pi0Boot=wl/n/(1-lambdas.get(i));
                double tempmse=mse.get(i);
                tempmse=tempmse+(pi0Boot-minPi0)*(pi0Boot-minPi0);
                mse.set(i,tempmse);
            }
        }
        double min=mse.get(0);
        int index=0;
        for (int i = 1; i < mse.size(); i++) {
            if(mse.get(i)<min){
                min=mse.get(i);
                index=i;
            }
        }
        double pi0=Math.max(Math.min(pi0s.get(index),1.0),0.0);
        return pi0;
    }


    /**
     * Method from Kall et al. used in q-value computation
     * @param list
     * @param maxSize
     * @return
     */
    private static ArrayList<Double> bootstrap(double[] list,int maxSize){
        ArrayList<Double> out=new ArrayList<>();
        int n=list.length;
        int numDraw=Math.min(n,maxSize);
        for (int i = 0; i < numDraw; i++) {
            int draw=(int)(rand.nextInt(n));
            //double pom=(double)PseudoRandom.nexRand()/((double)PseudoRandom.getRandMax()+1.0);
            //System.out.println("pom "+pom);
            //int draw=(int)(pom*n);
            out.add(list[draw]);
        }
        Collections.sort(out);
        //Collections.reverse(out);  //Not sure, they are claiming descending but the code looks like ascending
        return out;
    }


    /**
     * Method from Kall et al. used in q-value computation
     * @param array
     * @param start
     * @param value
     * @return
     */
    private static int lowerBound(double[] array,int start,double value){
        int pointer=array.length;
        for (int i = start; i < array.length; i++) {
            if(value<=array[i]){
                pointer=i;
                break;
            }
        }
        return pointer;
    }


    /**
     * Method from Kall et al. used in q-value computation
     * @param array
     * @param start
     * @param value
     * @return
     */
    private static int lowerBound(List<Double> array,int start,double value){
        int pointer=array.size();
        for (int i = start; i < array.size(); i++) {
            if(value>=array.get(i)){
                pointer=i;
                break;
            }
        }
        return pointer;
    }

    /**
     * Identify gene pair corresponding to index of p-value in a list of p-values.
     * @param x index of a p-value in a list of p-values created by methods in this class
     * @return indices of two gene such that the p-value corresponds to this pair
     */
    public int[] getIndicesToOrder(int x){
        int K=genesNumber*(genesNumber-1)/2;
        int k=(int)Math.ceil((Math.sqrt(8*(K-x))-1)/2.0);
        int[] coordinates=new int[2];  //row, column
        int difference=x-(K-(k*(k+1))/2);
        int rowNumber=genesNumber-1-k;
        int columnNumber=rowNumber+1+difference;
        coordinates[0]=rowNumber;
        coordinates[1]=columnNumber;
        return coordinates;
    }



}
