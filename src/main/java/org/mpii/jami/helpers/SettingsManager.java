package org.mpii.jami.helpers;

import java.util.HashMap;

/**
 * Created by mlist on 10/27/17.
 */
public class SettingsManager {

    private HashMap<String, Object> settings;
    public static final int defaultBatchSize = 100;
    public static final int defaultPermutations = 1000;
    public static final String defaultMethod = "";
    public static final int defaultThreads = -1;

    public SettingsManager(){
        this.settings = new HashMap<>();

        settings.put("header", true);
        settings.put("method", defaultMethod);
        settings.put("numberOfBins", 0);
        settings.put("tripleFormat", true);
        settings.put("numberOfThreads", defaultThreads);
        settings.put("numberOfPermutations", defaultPermutations);
        settings.put("pValueCutoff", 1.0);
        settings.put("selectedGenes", null);
        settings.put("restricted", false);
        settings.put("batchSize", defaultBatchSize);
        settings.put("considerZeros", true);
    }

    public SettingsManager(HashMap<String, Object> settings){
        this.settings = settings;
    }

    public Object get(String option){
        return(settings.get(option));
    }

    public void set(String option, Object obj){
        settings.put(option, obj);
    }
}
