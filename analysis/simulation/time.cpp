/***********************************************
 * time study from AHDC::hits
 *
 * on request of Michael
 *
 * @author Felix Touchte Codjo
 * @date October 7, 2025
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
int layer2number(int digit);
// end utilities

int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    

    //printf("argc : %d\n", argc);
    //if (argc < 2) { 
    //    return printf("Please, provide a filename...\n");
    //}
    //const char* filename = argv[1];
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/time/new/new-rec-all-23003.0-5.hipo";
    const char* filename = "/home/touchte-codjo/Desktop/amon/analysis/simulation/elastics_events.hipo";
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
    TH1D* H1_time = new TH1D("time", "time; time (ns); count", 100, -15, 400);
    std::vector<TH1D*> VecH1_time;
    VecH1_time.push_back(new TH1D("time_layer=11", "time_layer=11; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=21", "time_layer=21; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=22", "time_layer=22; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=31", "time_layer=31; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=32", "time_layer=32; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=41", "time_layer=41; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=42", "time_layer=42; time (ns); count", 100, -15, 400));
    VecH1_time.push_back(new TH1D("time_layer=51", "time_layer=51; time (ns); count", 100, -15, 400));
    TH2D* H2_t2d = new TH2D("t2d", "time2distance; time (ns); distance (mm)", 100, -15, 400, 100, -0.1, 5);


    long unsigned int nevents =0;
    //AhdcCCDB ccdb_new;
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
        event.getStructure(hitBank);
        event.getStructure(trackBank);
        event.getStructure(runConfigBank);
        event.getStructure(recEventBank);
        
        double startTime = recEventBank.getFloat("startTime", 0);
        H1_startTime->Fill(startTime);
        
        // t0
        for (int i = 0; i < hitBank.getRows(); i++) {
            int layer = hitBank.getByte("layer", i);
            int superlayer = hitBank.getByte("superlayer", i);
            //int component = hitBank.getShort("wire", i);
            double time = hitBank.getDouble("time", i);
            double doca = hitBank.getDouble("doca", i);
            double residual = hitBank.getDouble("residual", i);
            int digit = superlayer*10 + layer;
            VecH1_time[layer2number(digit)-1]->Fill(time);
            H1_time->Fill(time);
            // residual = hit.doca - kf.doca
            //H2_t2d->Fill(time, residual - doca);
            H2_t2d->Fill(time, doca - residual);
        }
    }
    printf("nevents    : %ld \n", nevents);
    

    TFile *f = new TFile("./output/time_calib_study.root", "RECREATE");
    H2_t2d->Write("time2distance");
    H1_startTime->Write("startTime");
    H1_time->Write("time");
    for (int i = 0; i < 8; i++) {
        VecH1_time[i]->Write(VecH1_time[i]->GetName());
    }
///////
    TCanvas* canvas1 = new TCanvas("c1_time","c1 time",1800, 990);
    TLegend* legend = new TLegend(0.5,0.6,0.9,0.9);
    int i_max = 0;
    int count_max = VecH1_time[0]->GetBinContent(VecH1_time[0]->GetMaximumBin());
    for (int i = 0; i < 8; i++) {
        VecH1_time[i]->SetLineColor(i+1);
        VecH1_time[i]->SetLineWidth(2);
        int count = VecH1_time[i]->GetBinContent(VecH1_time[i]->GetMaximumBin());
        if (count > count_max) {
            i_max = i;
            count_max = count;
        }
    }
    gStyle->SetOptStat(0);
    VecH1_time[i_max]->Draw();
    legend->AddEntry(VecH1_time[i_max], VecH1_time[i_max]->GetName());
    VecH1_time[i_max]->SetTitle("leadingEdgeTime - t0 - startTime");
    for (int i = 0; i < 8; i++) {
        if (i != i_max) {
            VecH1_time[i]->Draw("same");
            legend->AddEntry(VecH1_time[i], VecH1_time[i]->GetName());
        }
    }
    legend->Draw();
    canvas1->Write("time_combined"); 
///////
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

int layer2number(int digit) {
	if      (digit == 11) {
		return 1;
	} 
	else if (digit == 21) {
		return 2;
	} 
	else if (digit == 22) {
		return 3;
	} 
	else if (digit == 31) {
		return 4;
	} 
	else if (digit == 32) {
		return 5;
	} 
	else if (digit == 41) {
		return 6;
	} 
	else if (digit == 42) {
		return 7;
	} 
	else if (digit == 51) {
		return 8;
	} else {
		return -1; // not a layer
	}
}
