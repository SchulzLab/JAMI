package org.mpii.jami.helpers;

import java.text.DecimalFormat;

/**
 * Created by mlist on 10/27/17.
 */
public final class Progressbar {

    public static void updateProgress(final double progressPercentage) {
        final int width = 50; // progress bar width in chars

        String s = "[";
        int i = 0;
        for (; i <= (int) (progressPercentage * width); i++) {
            s += ".";
        }
        for (; i < width; i++) {
            s += " ";
        }
        s += "](" + (new DecimalFormat("#0.00")).format(progressPercentage * 100) + "%)";
        System.out.print(s + "\r");
    }
}
