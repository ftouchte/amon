#include "fOptions.h"
#include <cstdio>

fOptions::fOptions(std::vector<std::string> _opts) { opts = _opts;}

bool fOptions::IsOption(std::string str) {
    for (std::string opt: opts) {
        if (opt.compare(str) == 0) {
            return true;
        }
    }
    return false;
}

void fOptions::LoadOptions(int argc, const char * argv[]) {
    for (int i = 1; i < argc; i++) {
        if (IsOption(argv[i]) && i < argc-1) {
            mopts[argv[i]] = !IsOption(argv[i+1]) ? argv[i+1] : "";
        }
    }
}

std::string fOptions::GetValue(std::string str) {
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
        printf("|   %s = %s\n", opt.c_str(), mopts[opt].c_str());
    }
}


