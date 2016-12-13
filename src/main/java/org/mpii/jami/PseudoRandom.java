package org.mpii.jami;

/**
 * Created by fuksova on 12/16/15.
 * Intended to be used in ScoreStructure to agree with original q-values computation. Not used now
 */
public class PseudoRandom {
    private static long seed=5721888l;
    private static long randMax=4294967291l;

    public static long nexRand(){
        //double pokus=seed*279470273.0;
        //System.out.println("pokus "+pokus);
        long temp=seed*279470273l;
        System.out.println("seed: "+seed);
        seed=(temp)%randMax;
        return seed;
    }

    public static void setSeed(long s){
        seed=s;
    }

    public static long getRandMax(){
        return randMax;
    }


}
