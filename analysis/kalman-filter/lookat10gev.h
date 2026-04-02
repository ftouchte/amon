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
    TH2D* H2_track_corr_p_dEdx; ///< correlation p versus dEdx for the tracks

    // --- hit
    TH1D* H1_hit_residual; ///< hit residuals

    // --- electron
    TH1D* H1_electron_vz; ///< electron vertex
    TH1D* H1_electron_p; ///< electron momentum
    TH1D* H1_electron_theta; ///< electron theta
    TH1D* H1_electron_phi; /// electron phi

    /// Constructor
    Histograms() {
        // --- track
        H1_track_vz = new TH1D("track_vz", "track vertex; vz (cm); count", 100, -16, 16);
        H1_track_p = new TH1D("track_p", "track momentum; p (GeV); count", 100, 0, 1.0);
        H1_track_theta = new TH1D("track_theta", "track theta; #theta (deg); count", 100, 0, 180);
        H1_track_phi = new TH1D("track_phi", "track phi; #phi (deg); count", 100, 0, 360);
        H1_track_nhits = new TH1D("track_nhits", "Number of hits per track; nhits; count", 17, 0, 17);
        H2_track_corr_p_dEdx = new TH2D("track_corr_p_dEdx", "track : dE/dx versus p; p (GeV); dE/dx (MeV/mm)", 50, 0, 1.0, 50, 0, 180);

        // --- phi
        H1_hit_residual = new TH1D("hit_residual", "residual ; residual (mm); count", 100, -3, 3);

        // --- electron
        H1_electron_vz = new TH1D("electron_vz", "electron vertex; vz (cm); count", 100, -16, 16);
        H1_electron_p = new TH1D("electron_p", "electron momentum; p (GeV); count", 100, 0, 10.6);
        H1_electron_theta = new TH1D("electron_theta", "electron theta; #theta (deg); count", 100, 0, 40);
        H1_electron_phi = new TH1D("electron_phi", "electron phi; #phi (deg); count", 100, 0, 360);

    }

    /// Destructor
    ~Histograms() {
        // -- track
        delete H1_track_vz;
        delete H1_track_p;
        delete H1_track_theta;
        delete H1_track_phi;
        delete H1_track_nhits;
        delete H2_track_corr_p_dEdx;

        // --- hit
        delete H1_hit_residual;
        
        // --- electron
        delete H1_electron_vz;
        delete H1_electron_p;
        delete H1_electron_theta;
        delete H1_electron_phi;
    }

    /// @brief Save as ROOT file
    void  SaveIn(TDirectory* dir) {
        dir->cd();

        // --- track
        H1_track_vz->Write("track_vz");
        H1_track_p->Write("track_p");
        H1_track_theta->Write("track_theta");
        H1_track_phi->Write("track_phi");
        H1_track_nhits->Write("track_nhits");
        H2_track_corr_p_dEdx->Write("track_corr_p_dEdx");

        // --- hit
        H1_hit_residual->Write("hit_residual");

        // --- electron
        H1_electron_vz->Write("electron_vz");
        H1_electron_p->Write("electron_p");
        H1_electron_theta->Write("electron_theta");
        H1_electron_phi->Write("electron_phi");

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
