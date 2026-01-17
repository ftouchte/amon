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
    printf("| Test GetValue(\"%s\") = %s\n", "-i", OPT.GetValue("-i").c_str());
}