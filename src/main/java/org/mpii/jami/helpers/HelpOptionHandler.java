package org.mpii.jami.helpers;

import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.OptionDef;
import org.kohsuke.args4j.spi.BooleanOptionHandler;
import org.kohsuke.args4j.spi.Setter;

/**
 * Created by mlist on 11/9/17.
 */
public class HelpOptionHandler extends BooleanOptionHandler {
    public HelpOptionHandler(CmdLineParser parser, OptionDef option, Setter<Boolean> setter) {
        super(parser, option, setter);
    }
    public String printDefaultValue() {
        return null;  // this prevents the default value to be printed in usage info
    }
}