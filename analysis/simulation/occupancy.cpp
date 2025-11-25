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
    
    const char* output = "./output/occupancy_study.root";

    printf("argc : %d\n", argc);
    if (argc < 2) { 
        return printf("Please, provide a filename...\n");
    }
    const char* filename = argv[1];
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/raw-D2-run-23003.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::event event;
    long unsigned int nevents =0;
    // Histograms
    std::vector<TH1D*> VecH1_occupancy;
    VecH1_occupancy.push_back(new TH1D("occupancy_all", "occupancy_all; wire; count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_all&raw_cuts", "occupancy_all&raw_cuts; wire; count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_wfType<=2", "occupancy_wfType; wire; count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_wfType<=2&cuts", "occupancy_wfType<=2&cuts; wire; count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_wfType<=2&cuts-strong", "occupancy_wfType<=2&cuts-strong; wire; count [%]", 576, 0, 576));

    //signature // AhdcCCDB::AhdcCCDB(std::string _connection, int _runNo, std::string _variation, std::string _timestamp)
    AhdcCCDB ahdcConstants("mysql://clas12reader@clasdb.jlab.org/clas12", 23003, "default", "2025-04-18_13-44-11");
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
        
        // Loop over waveforms
        for (int i = 0; i < adcBank.getRows(); i++) {
            double sector    = adcBank.getInt("sector", i);
            double layer     = adcBank.getInt("layer", i);
            double component = adcBank.getInt("component", i);
            double leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", i);
            //double timeMax         = adcBank.getFloat("time", i);
            double timeOverThreshold = adcBank.getFloat("timeOverThreshold", i);
            int    wfType  = adcBank.getInt("wfType", i);
            int    adcMax  = adcBank.getInt("ADC", i);
            int    adcOffset  = adcBank.getFloat("ped", i);
            double t0         = ahdcConstants.get_t0(sector, layer, component).t0;
            double time       = leadingEdgeTime - t0;
            int nsamples = 0;
            for (int s = 0; s < 30; s++) {
                char buffer[50];
                sprintf(buffer, "s%d", s +1);
                int value = wfBank.getShort(buffer, i);
                if (value > 0) nsamples++;
            } 
            ahdcRawCuts rawCuts = ahdcConstants.get_rawCuts(sector, layer, component);
            double t_min = rawCuts.t_min;
            double t_max = rawCuts.t_max;
            double tot_min = rawCuts.tot_min;
            double tot_max = rawCuts.tot_max;
            double adc_min = rawCuts.adc_min;
            double adc_max = rawCuts.adc_max;
            double ped_min = rawCuts.ped_min;
            double ped_max = rawCuts.ped_max;
            // all
            VecH1_occupancy[0]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
            // rawCuts
            if ((adcMax >= adc_min) && (adcMax <= adc_max) && (adcOffset >= ped_min) && (adcOffset <= ped_max) && (leadingEdgeTime >= t_min) && (leadingEdgeTime <= t_max) && (timeOverThreshold >= tot_min) && (timeOverThreshold <= tot_max)) {
                VecH1_occupancy[1]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
            }
            // wfTye <= 2, then cuts
            if (wfType <= 2) {
                VecH1_occupancy[2]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
                if ((time >= 0) && (time <= 300) && (timeOverThreshold >= 300) && (timeOverThreshold <= 750) && (adcOffset >= 120) && (adcOffset <= 350) && (nsamples > 10) && (adcMax >= 0) && (adcMax <= 4095)) {
                    VecH1_occupancy[3]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
                }
                if ((time >= 0) && (time <= 300) && (timeOverThreshold >= 340) && (timeOverThreshold <= 620) && (adcOffset >= 120) && (adcOffset <= 350) && (nsamples > 14) && (adcMax >= 0) && (adcMax <= 4095)) {
                    VecH1_occupancy[4]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
                }
            }
        }
    }
    printf("nevents    : %ld \n", nevents);
    
    // Renormalise occupancy
    for (int i = 0; i < 5; i++) {
        VecH1_occupancy[i]->SetLineWidth(3);
        VecH1_occupancy[i]->Scale(100.0/nevents);
    }
    TCanvas* canvas1 = new TCanvas("c1_occupacy", "Occupancy study;wire; count [%]");
    // all
    VecH1_occupancy[0]->SetTitle("Occupancy");
    VecH1_occupancy[0]->SetLineColor(kBlack);
    VecH1_occupancy[0]->Draw();
    // raw cuts
    VecH1_occupancy[1]->SetLineColor(kBlue);
    VecH1_occupancy[1]->Draw("same");
    // wfType <= 2
    VecH1_occupancy[2]->SetLineColor(kGreen);
    VecH1_occupancy[2]->Draw("same");
    // wfType <=2 & raw cuts
    VecH1_occupancy[3]->SetLineColor(kRed);
    VecH1_occupancy[3]->Draw("same");
    // legend
    TLegend* legend = new TLegend(0.5,0.6,0.9,0.9);
    legend->AddEntry(VecH1_occupancy[0], VecH1_occupancy[0]->GetName());
    legend->AddEntry(VecH1_occupancy[2], VecH1_occupancy[2]->GetName());
    legend->AddEntry(VecH1_occupancy[3], VecH1_occupancy[3]->GetName());
    legend->AddEntry(VecH1_occupancy[1], VecH1_occupancy[1]->GetName());
    legend->Draw();

    // output
    TFile *f = new TFile(output, "RECREATE");
   
    std::vector<TDirectory*> wfType_DIR;

    // Loop over type
    canvas1->Write("COMBINED_occupancy");
    for (int i = 0; i < 5; i++) {
            VecH1_occupancy[i]->Write(VecH1_occupancy[i]->GetName());
    }

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

