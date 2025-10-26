/***********************************************
 * Study the Kalman Filter for the AHDC
 *
 * Distance in simulation
 *
 * @author Felix Touchte Codjo
 * @date October 24, 2025
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
    
    //const char * output = "./output/distance.root";
    const char * output = "./output/distance-v7.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v2.hipo";
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v7.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    //hipo::bank  trackBank(factory.getSchema("AHDC::track"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nMCtracks =0;
    long unsigned int nKFtracks =0;
    // Histograms
    // hit 
    TH1D* H1_mcMeanTime  = new TH1D("mcMeanTime", "mcMeanTime; mcMeanTime (ns); count", 100, 0, 250);
    TH1D* H1_decodedTime = new TH1D("decodedTime", "decodedTime; decodedTime (ns); count", 100, 0, 250);
    TH1D* H1_deltaTime = new TH1D("deltaTime", "deltaTime = mcMeanTime - decodedTime; deltaTime (ns); count", 100, -50, 50);
    TH1D* H1_mcMeanDistance  = new TH1D("mcMeanDistance", "mcMeanDistance; mcMeanDistance (mm); count", 100, 0, 4);
    TH1D* H1_decodedDistance = new TH1D("decodedDistance", "decodedDistance; decodedDistance (mm); count", 100, 0, 4);
    TH1D* H1_trueDistance    = new TH1D("trueDistance", "true distance (doca); trueDistance (mm); #count", 100, 0, 4); 
    std::vector<TH2D*> H2_corrDistances;
    H2_corrDistances.push_back(new TH2D("corr_1", "trueDistance vs decodedDistance;trueDistance (mm); decodedDistance (mm)", 100, 0, 4, 100, 0, 4));
    H2_corrDistances.push_back(new TH2D("corr_2", "trueDistance vs mcMeanDistance;trueDistance (mm); mcMeanDistance (mm)", 100, 0, 4, 100, 0, 4));
    H2_corrDistances.push_back(new TH2D("corr_3", "mcMeanDistance vs decodedDistance;mcMeanDistance (mm); decodedDistance (mm)", 100, 0, 4, 100, 0, 4));
    TH2D* H2_time2distance = new TH2D("corr_t2d", "mcMeanTime vs trueDistance;mcMeanTime (ns); trueDistance (mm)", 100, 0, 250, 100, 0, 4);
    TH1D* H1_t0  = new TH1D("t0", "t0; t0 (ns); count", 100, 0, 400);

    AhdcCCDB ahdcConstants;
    // Loop over events
    while( reader.next()){
        nevents++;
        // Progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(trackBank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);

        nMCtracks += mcBank.getRows();
        nKFtracks += trackBank.getRows();
        
        for (int i = 0; i < wfBank.getRows(); i++) {
            int sector = wfBank.getShort("sector", i);
            int layer = wfBank.getShort("layer", i);
            int component = wfBank.getShort("component", i);
            double mcMeanTime = 0.01*wfBank.getInt(std::string("s30").c_str(), i);
            double trueDistance = 0.001*wfBank.getInt(std::string("s28").c_str(), i);
            double t0 = ahdcConstants.get_t0(sector, layer, component).t0;
            double decodedTime = adcBank.getDouble("leadingEdgeTime", i) - t0; // for the simu, startTime = 0
            double mcMeanDistance = (mcMeanTime < 0) ? 0 : -0.0497 -0.00667*mcMeanTime + 0.389*sqrt(mcMeanTime) - 0.189*pow(mcMeanTime, 1.0/3);
            double decodedDistance = (decodedTime < 0) ? 0 : -0.0497 -0.00667*decodedTime + 0.389*sqrt(decodedTime) - 0.189*pow(decodedTime, 1.0/3);
            H1_mcMeanTime->Fill(mcMeanTime);
            H1_decodedTime->Fill(decodedTime);
            H1_deltaTime->Fill(mcMeanTime - decodedTime);
            H1_mcMeanDistance->Fill(mcMeanDistance);
            H1_decodedDistance->Fill(decodedDistance);
            H1_trueDistance->Fill(trueDistance);
            H1_t0->Fill(t0);
            H2_corrDistances[0]->Fill(trueDistance, decodedDistance);
            H2_corrDistances[1]->Fill(trueDistance, mcMeanDistance);
            H2_corrDistances[2]->Fill(mcMeanDistance, decodedDistance);
            H2_time2distance->Fill(mcMeanTime, trueDistance);
        }
    }
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    TFile *f = new TFile(output, "RECREATE");
    
    H1_mcMeanTime->Write("mcMeanTime");
    H1_decodedTime->Write("decodedTime");
    H1_mcMeanDistance->Write("mcMeanDistance");
    H1_decodedDistance->Write("decodedDistance");
    H1_trueDistance->Write("trueDistance");
    H2_corrDistances[0]->Write("corr_trueDistance_decodedDistance");
    H2_corrDistances[1]->Write("corr_trueDistance_mcMeanDistance");
    H2_corrDistances[2]->Write("corr_mcMeanDistance_decodedDistance");
    H2_time2distance->Write("time2distance");
    H1_t0->Write("t0");
    H1_deltaTime->Write("deltaTime");

    f->Close();
    printf("File created : %s\n", output);

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
