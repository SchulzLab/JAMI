package org.mpii.jami.model;

import java.util.function.Predicate;

/**
 * Created by mlist on 10/26/17.
 */
public class Triplet {

    private String geneOne;
    private String geneTwo;
    private String miRNA;
    private int geneOneIndex;
    private int geneTwoIndex;
    private int miRNAIndex;
    private double cmi;
    private double pValue;
    private double pAdjust;

    public Triplet(String[] args){
        this(args[0], args[1], args[2]);
    }

    public Triplet(String geneOne, String geneTwo, String miRNA) {
        this.geneOne = geneOne;
        this.geneTwo = geneTwo;
        this.miRNA = miRNA;
    }

    public int getGeneOneIndex() {
        return geneOneIndex;
    }

    public void setGeneOneIndex(int geneOneIndex) {
        this.geneOneIndex = geneOneIndex;
    }

    public int getGeneTwoIndex() {
        return geneTwoIndex;
    }

    public void setGeneTwoIndex(int geneTwoIndex) {
        this.geneTwoIndex = geneTwoIndex;
    }

    public int getMiRNAIndex() {
        return miRNAIndex;
    }

    public void setMiRNAIndex(int miRNAIndex) {
        this.miRNAIndex = miRNAIndex;
    }

    public double getCmi() {
        return cmi;
    }

    public void setCmi(double cmi) {
        this.cmi = cmi;
    }

    public double getpValue() {
        return pValue;
    }

    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    public String getGeneOne() {
        return geneOne;
    }

    public String getGeneTwo() {
        return geneTwo;
    }

    public String getMiRNA() {
        return miRNA;
    }

    public String[] toStringArray(){
        return new String[]{this.geneOne, this.geneTwo, this.miRNA};
    }

    private boolean isValid = true;

    public void markInvalid(){
        this.isValid = false;
    }

    @Override
    public int hashCode() {
        return geneOne.hashCode() + geneTwo.hashCode() + miRNA.hashCode();
    }

    public boolean equals(Object obj) {
        if(obj.getClass() != this.getClass()) return(false);
        Triplet triplet = (Triplet) obj;
        if(!this.geneOne.equals(triplet.geneOne)) return(false);
        if(!this.geneTwo.equals(triplet.geneTwo)) return(false);
        if(!this.miRNA.equals(triplet.miRNA)) return(false);
        return(true);
    }

    public static Predicate<Triplet> check(){
        return t -> t.isValid;
    };

    public void setpAdjust(Double pAdjust) {
        this.pAdjust = pAdjust;
    }

    public double getpAdjust() {
        return(this.pAdjust);
    }
}
