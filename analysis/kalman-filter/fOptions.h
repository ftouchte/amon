/***********************************************
 * Deal with options in C++
 *
 * @author Felix Touchte Codjo
 * @date January 17, 2026
 * ********************************************/

#ifndef F_OPTIONS_H
#define F_OPTIONS_H

#include <map>
#include <string>
#include <vector>

class fOptions {
    std::vector<std::string> opts;
    std::map<std::string, std::vector<std::string>> mopts;
public:
    /**
     * @brief Construct a new fOptions object
     * 
     * The user is free to give any sense to the list of options. 
     * e.g fOptions({"-i", "-o"}) where -i stands for input and -o for output
     * 
     * @param _opts list of options
     */
    fOptions(std::vector<std::string> _opts);

    /**
     * @brief Load/parse options
     * 
     * @param argc same argument as for main()
     * @param argv same argument as for main()
     */
    void LoadOptions(int argc, const char * argv[]);
    
    /**
     * @brief Check if str is an option
     * 
     * @param str 
     * @return true str is an option defined at the construction of the object
     * @return false otherwise
     */
     bool IsOption(std::string str);
    
     /**
     * @brief Get the first value of the option
     * 
     * @param str 
     * @return argument for this option
     */
    std::string GetValue(std::string str);

    /**
     * @brief Get all values for this option
     * 
     * @param str 
     * @return list of arguments for this option
     */
    std::vector<std::string> GetValues(std::string str);
    
    /**
     * @brief Display the value for each option
     * 
     */
    void Show();
};

#endif