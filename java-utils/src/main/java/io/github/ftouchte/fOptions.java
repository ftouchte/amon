package io.github.ftouchte;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/** Class for managing options in java
 * 
 * java adaptation of amon/analysis/kalman-filter/fOptions.cpp
 * 
 * @author Felix Touchte Codjo
 */
public class fOptions {
    ArrayList<String> opts = new ArrayList<>();
    Map<String, ArrayList<String>> mopts = new HashMap<>();

    /**
     * @brief Construct a new fOptions object
     * 
     * The user is free to give any sense to the list of options. 
     * e.g fOptions("-i", "-o") where -i stands for input and -o for output
     * 
     * @param _opts list of options
     */
    fOptions(String... _opts) {
        for (String x : _opts) {
            opts.add(x);
            mopts.put(x, new ArrayList<>());
        }
    }

    /**
     * @brief Load/parse options
     * 
     * @param args
     */
    void LoadOptions(String[] args) {
        String opt = "default";
        mopts.put(opt, new ArrayList<>());
        for (String x : args) {
            if (IsOption(x)) {
               opt = x;
            }
            else {
                mopts.get(opt).add(x);
            }
        }
    }

    /**
     * @brief Check if str is an option
     * 
     * @param str 
     * @return true str is an option defined at the construction of the object
     * @return false otherwise
     */
    boolean IsOption(String str) {
        for (String x : opts) {
            if (x.equals(str)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Get all values for this option
     * 
     * @param str 
     * @return list of arguments for this option
     */
    ArrayList<String> GetValues(String str) {
        if (IsOption(str)) {
            return new ArrayList<>(mopts.get(str));
        }
        else {
            return null;
        }
        
    }

    /**
     * @brief Get the first value of the option
     * 
     * @param str
     * @return first argument for this option
     */
    String GetValue(String str) {
        if (IsOption(str)) {
            if (mopts.get(str).size() > 0) {
                return mopts.get(str).get(0);
            }
            else {
                return "";
            }
        } else {
            return "";
        }
    }

    /**
     * @brief Display the value for each option
     * 
     */
    void Show() {
        int max_length = 0;
        for (String x : opts) {
            if (x.length() > max_length) {
                max_length = x.length();
            }
        }
        System.out.println("| List of options");
        for (String opt : opts) {
            System.out.printf("|   %-" + max_length + "s = ", opt);
            int counter = 0;
            ArrayList<String> list = GetValues(opt);
            for (String x : list) {
                if (counter >= 10) break;
                System.out.printf("%s ", x);
            }
            if (counter < 10) {
                System.out.println();
            } else {
                System.out.printf("and %d more...\n", list.size() - counter);
            }
        }
    }




}
