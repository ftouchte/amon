/***********************************************
 * Debug hit selection for 10.6 GeV run
 *
 * @author Felix Touchte Codjo
 * @date September 07, 2026
 * ********************************************/

#ifndef DEBUG_HIT_SELECTION_H
#define DEBUG_HIT_SELECTION_H

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"

void progressBar(int state, int bar_length = 100);
int layer2number(int digit);
int slc2wire(int sector, int layer, int component);

struct Histograms {

    TH1D* H1_any_hit_time;
    TH1D* H1_any_hit_tot;
    TH1D* H1_any_hit_ped;
    TH1D* H1_any_hit_amplitude;
    TH1I* H1_any_hit_wfType;
    TH1D* H1_any_hit_occupancy;

    TH1D* H1_selected_hit_time;
    TH1D* H1_selected_hit_tot;
    TH1D* H1_selected_hit_ped;
    TH1D* H1_selected_hit_amplitude;
    TH1I* H1_selected_hit_wfType;
    TH1D* H1_selected_hit_occupancy;

    TH1D* H1_track_nhits;
    TH1D* H1_track_hit_time;
    TH1D* H1_track_hit_tot;
    TH1D* H1_track_hit_ped;
    TH1D* H1_track_hit_amplitude;
    TH1I* H1_track_hit_wfType;
    TH1D* H1_track_hit_occupancy;

    Histograms() {
        
        H1_any_hit_time = new TH1D("any_hit_time", "any hit time; time (ns); count", 100, 0, 300);
        H1_any_hit_tot = new TH1D("any_hit_tot", "any hit tot; time over threshold (ns); count", 100, 100, 800);
        H1_any_hit_ped = new TH1D("any_hit_ped", "any hit pedestal; pedestal (ADC); count", 100, 0, 600);
        H1_any_hit_amplitude = new TH1D("any_hit_amplitude", "any hit amplitude; amplitude (ADC); count", 100, 0, 4000);
        H1_any_hit_wfType = new TH1I("any_hit_wfType", "any hit time; time (ns); count", 7, 0, 7);
        H1_any_hit_occupancy = new TH1D("any_hit_occupancy", "any hit time; time (ns); count", 576, 0, 576);

        H1_selected_hit_time = new TH1D("selected_hit_time", "selected hit time; time (ns); count", 100, 0, 300);
        H1_selected_hit_tot = new TH1D("selected_hit_tot", "selected hit tot; time over threshold (ns); count", 100, 100, 800);
        H1_selected_hit_ped = new TH1D("selected_hit_ped", "selected hit pedestal; pedestal (ADC); count", 100, 0, 600);
        H1_selected_hit_amplitude = new TH1D("selected_hit_amplitude", "selected hit amplitude; amplitude (ADC); count", 100, 0, 4000);
        H1_selected_hit_wfType = new TH1I("selected_hit_wfType", "selected hit time; time (ns); count", 7, 0, 7);
        H1_selected_hit_occupancy = new TH1D("selected_hit_occupancy", "selected hit time; time (ns); count", 576, 0, 576);

        H1_track_hit_time = new TH1D("track_hit_time", "track hit time; time (ns); count", 100, 0, 300);
        H1_track_hit_tot = new TH1D("track_hit_tot", "track hit tot; time over threshold (ns); count", 100, 100, 800);
        H1_track_hit_ped = new TH1D("track_hit_ped", "track hit pedestal; pedestal (ADC); count", 100, 0, 600);
        H1_track_hit_amplitude = new TH1D("track_hit_amplitude", "track hit amplitude; amplitude (ADC); count", 100, 0, 4000);
        H1_track_hit_wfType = new TH1I("track_hit_wfType", "track hit time; time (ns); count", 7, 0, 7);
        H1_track_hit_occupancy = new TH1D("track_hit_occupancy", "track hit time; time (ns); count", 576, 0, 576);
        H1_track_nhits = new TH1D("track_nhits", "track nhits", 12, 0, 12);

    }

    ~Histograms() {
        delete H1_any_hit_time;
        delete H1_any_hit_tot;
        delete H1_any_hit_ped;
        delete H1_any_hit_amplitude;
        delete H1_any_hit_wfType;
        delete H1_any_hit_occupancy;

        delete H1_selected_hit_time;
        delete H1_selected_hit_tot;
        delete H1_selected_hit_ped;
        delete H1_selected_hit_amplitude;
        delete H1_selected_hit_wfType;
        delete H1_selected_hit_occupancy;

        delete H1_track_nhits;
        delete H1_track_hit_time;
        delete H1_track_hit_tot;
        delete H1_track_hit_ped;
        delete H1_track_hit_amplitude;
        delete H1_track_hit_wfType;
        delete H1_track_hit_occupancy;
    }


    // useful
    // 96,116s/delete \(.*\);/\1->Write(\1->GetName());/gc
    void WriteIn(TDirectory* dir) {

        H1_any_hit_time->Write(H1_any_hit_time->GetName());
        H1_any_hit_tot->Write(H1_any_hit_tot->GetName());
        H1_any_hit_ped->Write(H1_any_hit_ped->GetName());
        H1_any_hit_amplitude->Write(H1_any_hit_amplitude->GetName());
        H1_any_hit_wfType->Write(H1_any_hit_wfType->GetName());
        H1_any_hit_occupancy->Write(H1_any_hit_occupancy->GetName());

        H1_selected_hit_time->Write(H1_selected_hit_time->GetName());
        H1_selected_hit_tot->Write(H1_selected_hit_tot->GetName());
        H1_selected_hit_ped->Write(H1_selected_hit_ped->GetName());
        H1_selected_hit_amplitude->Write(H1_selected_hit_amplitude->GetName());
        H1_selected_hit_wfType->Write(H1_selected_hit_wfType->GetName());
        H1_selected_hit_occupancy->Write(H1_selected_hit_occupancy->GetName());

        H1_track_nhits->Write(H1_track_nhits->GetName());
        H1_track_hit_time->Write(H1_track_hit_time->GetName());
        H1_track_hit_tot->Write(H1_track_hit_tot->GetName());
        H1_track_hit_ped->Write(H1_track_hit_ped->GetName());
        H1_track_hit_amplitude->Write(H1_track_hit_amplitude->GetName());
        H1_track_hit_wfType->Write(H1_track_hit_wfType->GetName());
        H1_track_hit_occupancy->Write(H1_track_hit_occupancy->GetName());
    }
};



#endif
