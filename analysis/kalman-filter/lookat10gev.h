/**********************************************
 * 
 * @author Felix Touchte Codjo
 * @date April 02, 2026
 **********************************************/

#ifndef LOOK_AT_10_GEV
#define LOOK_AT_10_GEV

#include <cstdio>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"

/**
 * @brief Object to store all the histograms used in the analysis
 * 
 * Dedicate the main code to event processing. When adding new histogram, remember to create them (contructor) and to delete them (destructor)
 * 
 */
struct Histograms {

    // --- monitoring
    
    // --- track
    TH1D* H1_track_vz; ///< track vertex
    TH1D* H1_track_p; ///< track momentum
    TH1D* H1_track_theta; ///< track theta
    TH1D* H1_track_phi; /// track phi
    TH1D* H1_track_nhits; /// number of hits per track
    TH1D* H1_track_chi2; /// track chi2
    TH2D* H2_track_corr_p_dEdx; ///< correlation p versus dEdx for the tracks

    // --- hit
    TH1D* H1_hit_residual; ///< hit residuals form tracks

    // --- electron
    TH1D* H1_electron_vz; ///< electron vertex
    TH1D* H1_electron_vz_nocut; ///< electron vertex no cut
    TH1D* H1_electron_p; ///< electron momentum
    TH1D* H1_electron_theta; ///< electron theta
    TH1D* H1_electron_phi; ///< electron phi

    TH1D* H1_delta_vz;

    // --- AHDC::adc
    TH1D* H1_hit_adc; 
    TH1D* H1_hit_time; 
    TH1D* H1_hit_leadingEdgeTime; 
    TH1D* H1_hit_tot;

    TH2D* H2_hit_adc_vs_time;
    TH2D* H2_hit_adc_vs_leadingEdgeTime;
    TH2D* H2_hit_adc_vs_tot;
    TH2D* H2_hit_time_vs_tot;
    TH2D* H2_hit_corrected_adc_vs_time;

    TH2D* H2_hit_slope_vs_adc_wfType0;
    TH2D* H2_hit_slope_vs_adc_allTypes;
    TH1D* H1_diff_adc_wfType0;



    /// Constructor
    Histograms() {
        // --- track
        H1_track_vz = new TH1D("track_vz", "track vertex; vz (cm); count", 100, -25, 25);
        H1_track_p = new TH1D("track_p", "track momentum; p (GeV); count", 100, 0, 1.0);
        H1_track_theta = new TH1D("track_theta", "track theta; #theta (deg); count", 100, 0, 180);
        H1_track_phi = new TH1D("track_phi", "track phi; #phi (deg); count", 100, 0, 360);
        H1_track_nhits = new TH1D("track_nhits", "Number of hits per track; nhits; count", 17, 0, 17);
        H1_track_chi2 = new TH1D("track_chi2", "chi2; chi2; count", 100, 0, 15);
        H2_track_corr_p_dEdx = new TH2D("track_corr_p_dEdx", "track : dE/dx versus p; p (GeV); dE/dx", 50, 0, 2.0, 50, 0, 500);

        // --- phi
        H1_hit_residual = new TH1D("hit_residual", "residual ; residual (mm); count", 100, -3, 3);

        // --- electron
        H1_electron_vz = new TH1D("electron_vz", "electron vertex; vz (cm); count", 100, -40, 25);
        H1_electron_vz_nocut = new TH1D("electron_vz_nocut", "electron vertex; vz (cm); count", 100, -60, 45);
        H1_electron_p = new TH1D("electron_p", "electron momentum; p (GeV); count", 100, 0, 10.6);
        H1_electron_theta = new TH1D("electron_theta", "electron theta; #theta (deg); count", 100, 0, 40);
        H1_electron_phi = new TH1D("electron_phi", "electron phi; #phi (deg); count", 100, 0, 360);

        H1_delta_vz = new TH1D("delta_vz", "delta vertex; #Delta vz (cm); count", 100, -16, 16);

        H1_hit_adc = new TH1D("hit_adc", "ADC; ADC; count", 100, 0, 4000); 
        H1_hit_time = new TH1D("hit_time", "time; time (ns); count", 100, 0, 300);
        H1_hit_leadingEdgeTime = new TH1D("hit_leadingEdgeTime", "leadingEdgeTime ; leadingEdgeTime (ns); count", 100, 200, 650);
        H1_hit_tot = new TH1D("hit_tot", "timeOverThreshold ; ToT (ns); count", 100, 200, 700);
        H2_hit_adc_vs_time = new TH2D("hit_adc_vs_time", "ADC vs time; ADC; time (ns)", 50, 0, 4000, 50, 0, 300);
        H2_hit_adc_vs_leadingEdgeTime = new TH2D("hit_adc_vs_leadingEdgeTime", "ADC vs leadingEdgeTime; ADC; leadingEdgeTime (ns)", 50, 0, 4000, 50, 200, 650);
        H2_hit_adc_vs_tot = new TH2D("hit_adc_vs_tot", "ADC vs timeOverThreshold; ADC; ToT (ns)", 50, 0, 4000, 50, 250, 700);
        H2_hit_corrected_adc_vs_time = new TH2D("hit_corrected_adc_vs_time", "ADC vs time including ToT correction; ADC; time (ns)", 50, 0, 4000, 50, 0, 300);
        H2_hit_time_vs_tot = new TH2D("hit_time_vs_tot", "time vs TimeOverthreshold; time (ns); toT (ns)", 50, 0, 300, 50, 200, 700);

        H2_hit_slope_vs_adc_wfType0 = new TH2D("hit_slope_vs_adc_wfType0", "slope vs ADC for wfType = 0; ADC; slope (ADC/ns)", 50, 0, 4000, 50, 0, 20);

        H2_hit_slope_vs_adc_allTypes = new TH2D("hit_slope_vs_adc_allTypes", "slope vs ADC all types; ADC; slope (ADC/ns)", 50, 0, 4000, 50, 0, 40);

        H1_diff_adc_wfType0 = new TH1D("hit_diff_adc_wfType0", "#Delta ADC = raw ADC minus predicted ADC;#Delta ADC; count", 100, -400, 400);

    }

    /// Destructor
    ~Histograms() {
        // -- track
        delete H1_track_vz;
        delete H1_track_p;
        delete H1_track_theta;
        delete H1_track_phi;
        delete H1_track_nhits;
        delete H1_track_chi2;
        delete H2_track_corr_p_dEdx;
        delete H1_delta_vz;
        delete H1_hit_residual;       
        
        // --- electron
        delete H1_electron_vz;
        delete H1_electron_vz_nocut;
        delete H1_electron_p;
        delete H1_electron_theta;
        delete H1_electron_phi;

        

        delete H1_hit_adc; 
        delete H1_hit_time; 
        delete H1_hit_leadingEdgeTime; 
        delete H1_hit_tot;
        delete H2_hit_adc_vs_time;
        delete H2_hit_adc_vs_leadingEdgeTime;
        delete H2_hit_adc_vs_tot;
        delete H2_hit_time_vs_tot;
        delete H2_hit_corrected_adc_vs_time;

        delete H2_hit_slope_vs_adc_wfType0;

        delete H2_hit_slope_vs_adc_allTypes;

        delete H1_diff_adc_wfType0;

        
    }

    /// @brief Save as ROOT file
    void  SaveIn(TDirectory* dir) {
        dir->cd();

        // delete \(.*\);/\1->Write(\1->GetName());/g

        // --- track
        TDirectory * track_dir = dir->mkdir("track");
        track_dir->cd();
        H1_track_vz->Write(H1_track_vz->GetName());
        H1_track_p->Write(H1_track_p->GetName());
        H1_track_theta->Write(H1_track_theta->GetName());
        H1_track_phi->Write(H1_track_phi->GetName());
        H1_track_nhits->Write(H1_track_nhits->GetName());
        H1_track_chi2->Write(H1_track_chi2->GetName());
        H2_track_corr_p_dEdx->Write(H2_track_corr_p_dEdx->GetName());
        H1_delta_vz->Write(H1_delta_vz->GetName());
        H1_hit_residual->Write(H1_hit_residual->GetName());   

        // --- electron
        TDirectory * electron_dir = dir->mkdir("electron");
        electron_dir->cd();
        H1_electron_vz->Write(H1_electron_vz->GetName());
        H1_electron_vz_nocut->Write(H1_electron_vz_nocut->GetName());
        H1_electron_p->Write(H1_electron_p->GetName());
        H1_electron_theta->Write(H1_electron_theta->GetName());
        H1_electron_phi->Write(H1_electron_phi->GetName());

        // --- AHDC::adc
        TDirectory * adc_dir = dir->mkdir("time_versus_adc");
        adc_dir->cd();

        H1_hit_adc->Write(H1_hit_adc->GetName()); 
        H1_hit_time->Write(H1_hit_time->GetName()); 
        H1_hit_leadingEdgeTime->Write(H1_hit_leadingEdgeTime->GetName()); 
        H1_hit_tot->Write(H1_hit_tot->GetName());
        H2_hit_adc_vs_time->Write(H2_hit_adc_vs_time->GetName());
        H2_hit_adc_vs_leadingEdgeTime->Write(H2_hit_adc_vs_leadingEdgeTime->GetName());
        H2_hit_adc_vs_tot->Write(H2_hit_adc_vs_tot->GetName());
        H2_hit_time_vs_tot->Write(H2_hit_time_vs_tot->GetName());
        H2_hit_corrected_adc_vs_time->Write(H2_hit_corrected_adc_vs_time->GetName());

        H2_hit_slope_vs_adc_wfType0->Write(H2_hit_slope_vs_adc_wfType0->GetName());

        H2_hit_slope_vs_adc_allTypes->Write(H2_hit_slope_vs_adc_allTypes->GetName());

        H1_diff_adc_wfType0->Write(H1_diff_adc_wfType0->GetName());
        

    }
};




struct futils {

    static void progressBar(int state, int bar_length = 100) { // state is a number between 0 and 100
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
};



#endif
