package org.mpii.jami.model;

/**
 * Created by mlist on 10/26/17.
 */
public class Triplet {

    private String geneOne;
    private String geneTwo;
    private String miRNA;
    private double cmi;
    private double pValue;

    public Triplet(String[] args){
        this(args[0], args[1], args[2]);
    }

    public Triplet(String geneOne, String geneTwo, String miRNA) {
        this.geneOne = geneOne;
        this.geneTwo = geneTwo;
        this.miRNA = miRNA;
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
}
