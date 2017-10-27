package org.mpii.jami.scores;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;

/**
 * Created by fuksova on 1/14/16.
 * This class contains some functions to show how to work with Cupid score file and how to use gene pair selection
 * as they use it. Those methods are not suitable for multiple reasons:
 * 1) p-values for evaluation gene pair significance are obtained from heuristic contingency tables on which Fisher
 * exact test is applied.
 * 2) Pairs for which most interactions scores are zero are marked as the most significant.
 * 3) We do not know whether the score values are reliable. According to several other incorrect heuristics
 * or even wrong usages of methods, we can assume that the score is not reliable.
 *
 */
public class GenePairsSelection {


    /**
     * How to work with score file, obtain p-values, select them according to threshold.
     */
    public static void prunnedPValuesComplete(){
        ScoreStructure.geneInd=1;  //index of column with gene name
        ScoreStructure.miRNAInd=2;  //index of column with miRNA name
        ScoreStructure.valueInd=0;  //index of column with score value
        //Reading modified score files of Cupid, those files contain only columns with names and score values
        ScoreStructure scoreStructure=new ScoreStructure("tablePosColumns.txt","tableNegColumns.txt");
        scoreStructure.init(true,true);
        double pValue=scoreStructure.computePValue("LRP12","MFF");  //How to compute p-value for one pair
        double [] pvalues=scoreStructure.getAllPValues();   //Computation of p-values from score values
        scoreStructure.pairsBellowThresholdToFile("pairsBelowThreshold2.txt",pvalues,0.00008706713543752097);  //selection of pairs below threshold and storing them into a file
    }


    /**
     * For debugging purposes
     */
    public static void inspectContingencyTables(){
        ScoreStructure.geneInd=1;
        ScoreStructure.miRNAInd=2;
        ScoreStructure.valueInd=0;
        ScoreStructure scoreStructure=new ScoreStructure("tablePosColumns.txt","tableNegColumns.txt");
        scoreStructure.init(true,true);
        double[][] tOrig=scoreStructure.computePseudoContTable(scoreStructure.geneNameToInt.get("ESR1"),scoreStructure.geneNameToInt.get("HIF1A"));
        LinkedList<String[]> pairs=new LinkedList<>();
        try {
            BufferedReader br=new BufferedReader(new FileReader("tightPairsBelowThreshold.txt"));
            br.readLine();
            String line=br.readLine();

            while (line!=null){
                String[] spl=line.split("\t");
                pairs.add(spl);
                line=br.readLine();

            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        LinkedList<double[][]> tables=new LinkedList<>();
        for (String[] pair : pairs) {
            double[][] t1=scoreStructure.computePseudoContTable(scoreStructure.geneNameToInt.get(pair[0]),scoreStructure.geneNameToInt.get(pair[1]));
            tables.add(t1);
        }
        System.out.println("done");
    }


    /**
     * For debugging purposes. Which pairs have the lowest p-value?
     */
    public static void writeGenePairsBellowP(){
        ScoreStructure.geneInd=1;
        ScoreStructure.miRNAInd=2;
        ScoreStructure.valueInd=0;
        ScoreStructure scoreStructure=new ScoreStructure("tablePosColumns.txt","tableNegColumns.txt");
        scoreStructure.init(true,true);
        scoreStructure.scoreMatrix=null;
        double experimentalyFoundPValue=0.00008706713543752097;
        double myTightPvalue=Math.pow(10,-29);
        double[] pv = ScoreStructure.getPValuesFromFile("pvalues.txt", 163669278);
        scoreStructure.pairsBellowThresholdToFile("tightPairsBelowThreshold.txt",pv,myTightPvalue);
    }


    /**
     * Copy of first n lines of file p-values
     */
    public static void rewritePValues(){
        double[] pv = ScoreStructure.getPValuesFromFile("pvalues.txt", 3531);
        ScoreStructure.writePValuesToFile("pvalues2.txt", pv);
    }

    /**
     * Some examples how to work with p-values, obtain q-values, and p-value threshold from q-values
     */
    public static void qValues(){
        //  PseudoRandom.setSeed(1457212008l);
        ScoreStructure.rand.setSeed(21545422l);
        double[] pv = ScoreStructure.getPValuesFromFile("pvalues.txt", 163669278);
        double[] qv = ScoreStructure.getQvaluesFromPvalues(pv,true);
        double pThreshold=ScoreStructure.getPThresholdToQThreshold(0.01, pv, "PQList.txt"); //This also contains computing of q-values from p-values
        System.out.println("p-value threshold: "+pThreshold);



    }

}
