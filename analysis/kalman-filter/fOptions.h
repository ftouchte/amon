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
    std::map<std::string, std::string> mopts;
public:
    fOptions(std::vector<std::string> _opts);
    void LoadOptions(int argc, const char * argv[]);
    bool IsOption(std::string str);
    std::string GetValue(std::string str);
    void Show();
    void RunTest(int argc, const char * argv[]);
};

#endif