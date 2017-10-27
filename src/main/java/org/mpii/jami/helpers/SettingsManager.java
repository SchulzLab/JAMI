package org.mpii.jami.helpers;

import java.util.HashMap;

/**
 * Created by mlist on 10/27/17.
 */
public class SettingsManager {

    private HashMap<String, Object> settings;

    public SettingsManager(){
        this.settings = new HashMap<>();

        settings.put("header", true);
        settings.put("method", "");
        settings.put("numberOfBins", 0);
        settings.put("tripleFormat", true);
        settings.put("numberOfThreads", -1);
        settings.put("numberOfPermutations", 1000);
        settings.put("pValueCutoff", 1.0);
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
