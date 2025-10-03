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
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/pull-request/extracted_new_simu.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/23003/dev/rec-run_23003.00000.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/all-rec-23003.hipo";
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
    hipo::bank  recEventBank(factory.getSchema("REC::Event"));
    hipo::event event;

    TH1D* H1_startTime = new TH1D("startTime", "startTime; startTime (ns); count", 100, 0, 170);
    TH1D* H1_t0_old = new TH1D("t0_old", "t0; t0 (ns); count", 100, 100, 400);
    TH1D* H1_t0_new = new TH1D("t0_new", "t0; t0 (ns); count", 100, 100, 400);
    TH1D* H1_time = new TH1D("time", "time; time (ns); count", 100, 0, 300);
    TH1D* H1_mctime = new TH1D("mctime", "mctime; mctime (ns); count", 100, 0, 300);
    TH1D* H1_doca = new TH1D("doca", "doca; doca (mm); count", 100, 0, 4);
    TH1D* H1_mcdoca = new TH1D("mcdoca", "mcdoca; mcdoca (mm); count", 100, 0, 4);
    TH1D* H1_diffdoca = new TH1D("diffdoca", "doca - mcdoca; diffdoca (mm); count", 100, -4, 4);
    TH1D* H1_difftime = new TH1D("difftime", "time - mctime; difftime (ns); count", 100, -50, 50);
    TH2D* H2_difftime_t0 = new TH2D("difftime_t0", "Delta time vs t0; t0 (ns); difftime (ns)", 100, 100, 250, 100, -50, 50);

    long unsigned int nevents =0;
    //AhdcCCDB ccdb_new;
    //AhdcCCDB ccdb_new("mysql://clas12reader@clasdb.jlab.org/clas12", 22092, "default", "2025-10-01_00-00-00");
    AhdcCCDB ccdb_new("mysql://clas12reader@clasdb.jlab.org/clas12", 22092, "default", "2025-10-02_00-00-00");
    //AhdcCCDB ccdb_old("mysql://clas12reader@clasdb.jlab.org/clas12", 23003, "default", "2025-09-28_00-00-00");
    AhdcCCDB ccdb_old("mysql://clas12reader@clasdb.jlab.org/clas12", 22092, "default", "2025-09-28_00-00-00");
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
        event.getStructure(recEventBank);
        
        //int evt = runConfigBank.getInt("event", 0);
        /*if (evt == 5272) {
            //adcBank.show();
            hitBank.show();
            trackBank.show();
            runConfigBank.show();
        }*/

        // startTime
        //double startTime = recEventBank.getFloat("startTime", 0);
        double startTime = 0; //mean value
        H1_startTime->Fill(startTime);
        
        // t0
        for (int i = 0; i < adcBank.getRows(); i++) {
            int sector = adcBank.getByte("sector", i);
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            double t0_old = ccdb_old.get_t0(sector, layer, component).t0;
            double t0_new = ccdb_new.get_t0(sector, layer, component).t0;
            H1_t0_old->Fill(t0_old);
            H1_t0_new->Fill(t0_new);
        
            // mc and rec distance
            double time = adcBank.getFloat("leadingEdgeTime", i) - t0_new - startTime;
            double mctime = wfBank.getShort("s30", i)*0.01;
            double doca = -0.0497 - 0.00667*pow(time,1.0) + 0.389*pow(time,1.0/2) - 0.189*pow(time,1.0/3);
            double mcdoca = -0.0497 - 0.00667*pow(mctime,1.0) + 0.389*pow(mctime,1.0/2) - 0.189*pow(mctime,1.0/3);
            H1_time->Fill(time);
            H1_mctime->Fill(mctime);
            H1_doca->Fill(doca);
            H1_mcdoca->Fill(mcdoca);
            H1_difftime->Fill(time-mctime);
            H1_diffdoca->Fill(doca-mcdoca);
            H2_difftime_t0->Fill(t0_new, time-mctime);
        }
        
    }
    printf("nevents    : %ld \n", nevents);
   
    TFile *f = new TFile("./output/fast_study.root", "RECREATE");
    H1_startTime->Write("startTime");
    H1_t0_old->Write("t0_old");
    H1_t0_new->Write("t0_new");
    H1_time->Write("time");
    H1_mctime->Write("mctime");
    H1_doca->Write("doca");
    H1_mcdoca->Write("mcdoca");
    H1_difftime->Write("difftime");
    H1_diffdoca->Write("diffdoca");
    H2_difftime_t0->Write("corr_difftime_t0");
    f->Close();

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

