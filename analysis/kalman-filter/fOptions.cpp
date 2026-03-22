/***********************************************
 * Deal with options in C++
 *
 * @author Felix Touchte Codjo
 * @date January 17, 2026
 * ********************************************/

#include "fOptions.h"
#include <cstdio>

fOptions::fOptions(std::vector<std::string> _opts) { 
    opts = _opts;
    for (auto x : _opts) {
        mopts[x] = {};
    }
}

bool fOptions::IsOption(std::string str) {
    for (std::string opt: opts) {
        if (opt.compare(str) == 0) {
            return true;
        }
    }
    return false;
}

void fOptions::LoadOptions(int argc, const char * argv[]) {
    // for (int i = 1; i < argc; i++) {
    //     if (IsOption(argv[i]) && i < argc-1) {
    //         mopts[argv[i]] = !IsOption(argv[i+1]) ? argv[i+1] : "";
    //     }
    // }
    std::string current_opt = "default";
    for (int i = 1; i < argc; i++) {
        if (IsOption(argv[i])) {
            current_opt = argv[i];
        } else {
            mopts[current_opt].push_back(argv[i]);
        }
    }
}

std::string fOptions::GetValue(std::string str) {
    if (mopts[str].size() > 0) {
        return mopts[str][0];
    }
    return "";
}

std::vector<std::string> fOptions::GetValues(std::string str) {
    return mopts[str];
}

void fOptions::Show() {
    int max_length = 0;
    for (std::string opt: opts) {
        if ((int) opt.size() > max_length) { 
            max_length = opt.size();
        }
    }
    printf("| List of options\n");
    for (std::string opt: opts) {
        //printf("|   %s = %s\n", opt.c_str(), mopts[opt].c_str());
        printf("|   %s = ", opt.c_str());
        for (auto x : mopts[opt]) {
            printf("%s ", x.c_str());
        }
        printf("\n");
    }
}


