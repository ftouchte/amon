/***********************************************
 * Study occupancy -- wfType cuts vs raw cuts
 *
 * @author Felix Touchte Codjo
 * @date September 26, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

#include "reader.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "Math/PdfFuncMathCore.h"
#include "TApplication.h"
#include "TText.h"
#include "THStack.h"

#include "AhdcCCDB.h"
#include "futils.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    

    //printf("argc : %d\n", argc);
    //if (argc < 2) { 
    //    return printf("Please, provide a filename...\n");
    //}
    //const char* filename = argv[1];
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/all-rec-23003.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/coat-13.0.1/rec_clas_023003.evio.00000.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/raw-D2-run-23003.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  runConfigBank(factory.getSchema("RUN::config"));
    hipo::event event;
    long unsigned int nevents =0;
    // Loop over events
    while( reader.next()){
        nevents++;
        // Progress Bar
        //if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
        //    progressBar(100.0*nevents/reader.getEntries());
        //}
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(hitBank);
        event.getStructure(trackBank);
        event.getStructure(runConfigBank);
        
        int evt = runConfigBank.getInt("event", 0);
        if (evt == 5272) {
            //adcBank.show();
            hitBank.show();
            trackBank.show();
            runConfigBank.show();
        }
        
    }
    printf("nevents    : %ld \n", nevents);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("* time elapsed : %lf seconds\n", elapsed.count());
    return 0;
}

void progressBar(int state, int bar_length) { // state is a number between 0 and 100
    // for the moment the bar length is not variable
    if (state > bar_length) {return ;}
    printf("\rProgress \033[32m\[");
    for (int i = 0; i <= state; i++) {
        printf("#");
    }
    printf("\033[0m");
    for (int i = state+1; i < bar_length; i++) {
        printf(".");
    }
    if (state == 100) {
        printf("\033[32m] \033[1m %d %%\033[0m\n", state);
    } else {
        printf("] %d %%", state);
    }
    fflush(stdout);
}

