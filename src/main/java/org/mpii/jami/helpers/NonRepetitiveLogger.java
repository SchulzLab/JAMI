package org.mpii.jami.helpers;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by mlist on 11/3/17.
 */
public class NonRepetitiveLogger {

    private static final Logger logger = LogManager.getLogger("JAMI");

    private Set<String> logMessages = new HashSet<>();

    public void warn(String message){
        if(!logMessages.contains(message)){
            logger.warn(message);
            logMessages.add(message);
        }
    }
}
