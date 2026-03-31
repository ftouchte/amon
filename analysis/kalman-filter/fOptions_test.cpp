#include "fOptions.h"
#include <cstdio>
#include <cassert>

int main(int argc, const char * argv[]) {
    printf("Test class : fOptions\n");
    fOptions OPT({"-i", "-o", "-v"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();
    assert(OPT.IsOption("-i"));
    assert(OPT.IsOption("-o"));
    assert(OPT.IsOption("-v"));
    // look at -i
    std::vector<std::string> vec = OPT.GetValues("-i");
    printf("| Test GetValue(\"-i\") = ");
    for (auto x : vec) {
        printf("%s ", x.c_str());
    }
    printf("\n");
}