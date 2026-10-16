/***********************************************
 * Study elastics events off He4
 *
 * @author Felix Touchte Codjo
 * @date October 15, 2026
 * 
 * ********************************************/



#ifndef HELIUM_ELASTICS_H
#define HELIUM_ELASTICS_H


#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"

using FourVector = ROOT::Math::PxPyPzEVector;
using ThreeVector = ROOT::Math::XYZVector;



struct particle
{
    FourVector momentum; ///< 4 momentum vector
    int row = -1; ///< CLAS bank row
};

struct kinematics {
    double Q2; ///< Q2 virtuality
    double xB; ///< x Bjorken
    double t; ///< momentum transferred to the hadron
    double Phi; ///< angle between the hadronic plane and the leptonic plane
    double W2;  ///< Missing mass of eX
    double y; ///< (E-E')/E
    double nu; ///< E-E'
};



struct Histograms
{

    TH1D* h1_Q2;
    TH1D* h1_W2;
    TH1D* h1_xB;
    TH1D* h1_nu;
    TH1D* h1_y;
    TH1D* h1_electron_vz_nocut;

    // after W2 cut and vz cut

    TH1D* h1_Q2_W2cut;
    TH1D* h1_W2_W2cut;
    TH1D* h1_xB_W2cut;
    TH1D* h1_nu_W2cut;
    TH1D* h1_y_W2cut;

    TH1D* h1_electron_p;
    TH1D* h1_electron_theta;
    TH1D* h1_electron_phi;
    TH1D* h1_electron_vz;

    // After selecting one electron and requiring good tracks (chi2 and nb_hits cut)
    TH1D* h1_ahdc_track_nhits_nocut;
    TH1D* h1_ahdc_track_chi2_nocut;
    TH1D* h1_ahdc_track_nhits;
    TH1D* h1_ahdc_track_chi2;
    TH1D* h1_ahdc_track_p;
    TH1D* h1_ahdc_track_theta;
    TH1D* h1_ahdc_track_phi;
    TH1D* h1_ahdc_track_vz;
    TH1D* h1_ahdc_track_delta_phi;

    // elastics
    TH1D* h1_elastics_ahdc_track_delta_phi;
    TH1D* h1_elastics_ahdc_track_nhits;
    TH1D* h1_elastics_ahdc_track_chi2;
    TH1D* h1_elastics_ahdc_track_p;
    TH1D* h1_elastics_ahdc_track_theta;
    TH1D* h1_elastics_ahdc_track_phi;
    TH1D* h1_elastics_ahdc_track_vz;
    TH2D* h2_elastics_ahdc_track_p_vs_dEdx;
    TH1D* h1_elastics_electron_p;
    TH1D* h1_elastics_electron_theta;
    TH1D* h1_elastics_electron_phi;
    TH1D* h1_elastics_electron_vz;
    TH1D* h1_elastics_delta_vz;
    TH1D* h1_elastics_expected_track_p;
    TH1D* h1_elastics_expected_track_theta;
    TH1D* h1_elastics_expected_track_phi;
    TH1D* h1_elastics_diff_p;
    TH1D* h1_elastics_diff_theta;
    TH1D* h1_elastics_diff_phi;
    

    /**
     * @brief Construct a new Histograms object
     * 
     */
    Histograms() {

        h1_Q2 = new TH1D("Q2", "Q2; Q2 (GeV^{2}); count", 100, 0, 0.5);
        h1_W2 = new TH1D("W2", "W2; W2 (GeV^{2}); count", 100, 13, 21);
        h1_xB = new TH1D("xB", "xB; xB; count", 100, 0, 1);
        h1_nu = new TH1D("nu", "#nu = E-E'; #nu (GeV); count", 100, 0, 1);
        h1_y =  new TH1D("y", "y = (E-E')/E; y; count", 100, 0, 0.4);
        h1_electron_vz_nocut = new TH1D("electron_vz_nocut", "electron vz; vz (cm); count", 100, -45, 30);

        h1_Q2_W2cut = new TH1D("Q2_W2cut", "Q2; Q2 (GeV^{2}); count", 100, 0, 0.2);
        h1_W2_W2cut = new TH1D("W2_W2cut", "W2; W2 (GeV^{2}); count", 100, 13.79, 14.05);
        h1_xB_W2cut = new TH1D("xB_W2cut", "xB; xB; count", 100, 0, 1);
        h1_nu_W2cut = new TH1D("nu_W2cut", "#nu = E-E'; #nu (GeV); count", 100, 0, 0.04);
        h1_y_W2cut =  new TH1D(" y_W2cut", "y = (E-E')/E; y; count", 100, 0, 0.05);

        h1_electron_p = new TH1D("electron_p", "electron p; p (GeV); count", 100, 2.19, 2.3);
        h1_electron_theta = new TH1D("electron_theta", "electron theta;  (deg); count", 100, 4, 14);
        h1_electron_phi = new TH1D("electron_phi", "electron phi; phi (deg); count", 100, 0, 360);
        h1_electron_vz = new TH1D("electron_vz", "electron vz; vz (cm); count", 100, -45, 30);

        h1_ahdc_track_nhits_nocut = new TH1D("ahdc_track_nb_hits_nocut", "ahdc_track nb. hits; nb. hits; count", 12, 0, 12);
        h1_ahdc_track_chi2_nocut = new TH1D("ahdc_track_chi2_nocut", "ahdc_track chi2; chi2; count", 100, 0, 10);

        h1_ahdc_track_nhits = new TH1D("ahdc_track_nb_hits", "ahdc_track nb. hits; nb. hits; count", 9, 6, 15);
        h1_ahdc_track_chi2 = new TH1D("ahdc_track_chi2", "ahdc_track chi2; chi2; count", 100, 0, 10);
        h1_ahdc_track_p = new TH1D("ahdc_track_p", "ahdc_track p; p (GeV); count", 100, 0, 1);
        h1_ahdc_track_theta = new TH1D("ahdc_track_theta", "ahdc_track theta;  (deg); count", 100, 0, 180);
        h1_ahdc_track_phi = new TH1D("ahdc_track_phi", "ahdc_track phi; phi (deg); count", 100, 0, 360);
        h1_ahdc_track_vz = new TH1D("ahdc_track_vz", "ahdc_track vz; vz (cm); count", 100, -25, 20);
        h1_ahdc_track_delta_phi = new TH1D("ahdc_track_delta_phi", "ahdc_track #Delta #Phi; #Delta #Phi (deg); count", 100, -30, 30);

        h1_elastics_ahdc_track_delta_phi = new TH1D("elastics_ahdc_track_delta_phi", "ahdc_track #Delta #Phi; #Delta #Phi (deg); count", 100, -8, 8);
        h1_elastics_ahdc_track_nhits = new TH1D("elastics_ahdc_track_nb_hits", "ahdc_track nb. hits; nb. hits; count", 9, 6, 15);
        h1_elastics_ahdc_track_chi2 = new TH1D("elastics_ahdc_track_chi2", "ahdc_track chi2; chi2; count", 100, 0, 10);
        h1_elastics_ahdc_track_p = new TH1D("elastics_ahdc_track_p", "ahdc_track p; p (GeV); count", 100, 0, 1);
        h1_elastics_ahdc_track_theta = new TH1D("elastics_ahdc_track_theta", "ahdc_track theta;  (deg); count", 100, 0, 180);
        h1_elastics_ahdc_track_phi = new TH1D("elastics_ahdc_track_phi", "ahdc_track phi; phi (deg); count", 100, 0, 360);
        h1_elastics_ahdc_track_vz = new TH1D("elastics_ahdc_track_vz", "ahdc_track vz; vz (cm); count", 100, -25, 20);
        h2_elastics_ahdc_track_p_vs_dEdx = new TH2D("elastics_ahdc_track_p_vs_dEdx", "ahdc track dEdx vs p", 50, 0, 1, 50, 0, 1500);
        h1_elastics_electron_p = new TH1D("elastics_electron_p", "electron p; p (GeV); count", 100, 2.19, 2.25);
        h1_elastics_electron_theta = new TH1D("elastics_electron_theta", "electron theta;  (deg); count", 100, 4, 14);
        h1_elastics_electron_phi = new TH1D("elastics_electron_phi", "electron phi; phi (deg); count", 100, 0, 360);
        h1_elastics_electron_vz = new TH1D("elastics_electron_vz", "electron vz; vz (cm); count", 100, -45, 30);
        h1_elastics_delta_vz = new TH1D("elastics_ahdc_track_delta_vz", "#Delta vz = vz^{electron} - vz^{AHDC}; #Delta vz (cm); count", 100, -11, 3);
        h1_elastics_expected_track_p = new TH1D("elastics_expected_track_p", "expected track p; p (GeV); count", 100, 0, 1);
        h1_elastics_expected_track_theta = new TH1D("elastics_expected_track_theta", "expected track theta;  (deg); count", 100, 77, 88);
        h1_elastics_expected_track_phi = new TH1D("elastics_expected_track_phi", "expected track phi; phi (deg); count", 100, 0, 360);

        h1_elastics_diff_p = new TH1D("elastics_diff_p", "elastics diff p; #Delta p (GeV); count", 100, -0.11, 0.11);
        h1_elastics_diff_theta = new TH1D("elastics_diff_theta", "elastics diff theta; #Delta #theta (deg); count", 100, -60, 60);
        h1_elastics_diff_phi = new TH1D("elastics_diff_phi", "elastics diff phi;  #Delta #phi (deg); count", 100, -8, 8);
        


    }

    /**
     * @brief Destroy the Histograms object
     * 
     */
    ~Histograms() {
        // 69,83s/TH[12]D\*/delete/g

        delete h1_Q2;
        delete h1_W2;
        delete h1_xB;
        delete h1_nu;
        delete h1_y;
        delete h1_electron_vz_nocut;

        delete h1_Q2_W2cut;
        delete h1_W2_W2cut;
        delete h1_xB_W2cut;
        delete h1_nu_W2cut;
        delete h1_y_W2cut;

        delete h1_electron_p;
        delete h1_electron_theta;
        delete h1_electron_phi;
        delete h1_electron_vz;

        delete h1_ahdc_track_nhits_nocut;
        delete h1_ahdc_track_chi2_nocut;
        delete h1_ahdc_track_nhits;
        delete h1_ahdc_track_chi2;
        delete h1_ahdc_track_p;
        delete h1_ahdc_track_theta;
        delete h1_ahdc_track_phi;
        delete h1_ahdc_track_vz;

        delete h1_ahdc_track_delta_phi;
           
        delete h1_elastics_ahdc_track_delta_phi;
        delete h1_elastics_ahdc_track_nhits;
        delete h1_elastics_ahdc_track_chi2;
        delete h1_elastics_ahdc_track_p;
        delete h1_elastics_ahdc_track_theta;
        delete h1_elastics_ahdc_track_phi;
        delete h1_elastics_ahdc_track_vz;
        delete h2_elastics_ahdc_track_p_vs_dEdx;
        delete h1_elastics_electron_p;
        delete h1_elastics_electron_theta;
        delete h1_elastics_electron_phi;
        delete h1_elastics_electron_vz;
        delete h1_elastics_delta_vz;
        delete h1_elastics_expected_track_p;
        delete h1_elastics_expected_track_theta;
        delete h1_elastics_expected_track_phi;

        delete h1_elastics_diff_p;
        delete h1_elastics_diff_theta;
        delete h1_elastics_diff_phi;
        

    }


    void Write(TDirectory* dir) {
        dir->cd();

        // 89,103s/TH[12]D\* h\(.*\);/h\1->Write(h\1->GetName());/g
        h1_Q2->Write(h1_Q2->GetName());
        h1_W2->Write(h1_W2->GetName());
        h1_xB->Write(h1_xB->GetName());
        h1_nu->Write(h1_nu->GetName());
        h1_y->Write(h1_y->GetName());
        h1_electron_vz_nocut->Write(h1_electron_vz_nocut->GetName());

        h1_Q2_W2cut->Write(h1_Q2_W2cut->GetName());
        h1_W2_W2cut->Write(h1_W2_W2cut->GetName());
        h1_xB_W2cut->Write(h1_xB_W2cut->GetName());
        h1_nu_W2cut->Write(h1_nu_W2cut->GetName());
        h1_y_W2cut->Write(h1_y_W2cut->GetName());

        h1_electron_p->Write(h1_electron_p->GetName());
        h1_electron_theta->Write(h1_electron_theta->GetName());
        h1_electron_phi->Write(h1_electron_phi->GetName());
        h1_electron_vz->Write(h1_electron_vz->GetName());

        h1_ahdc_track_nhits_nocut->Write(h1_ahdc_track_nhits_nocut->GetName());
        h1_ahdc_track_chi2_nocut->Write(h1_ahdc_track_chi2_nocut->GetName());
        h1_ahdc_track_nhits->Write(h1_ahdc_track_nhits->GetName());
        h1_ahdc_track_chi2->Write(h1_ahdc_track_chi2->GetName());

        h1_ahdc_track_p->Write(h1_ahdc_track_p->GetName());
        h1_ahdc_track_theta->Write(h1_ahdc_track_theta->GetName());
        h1_ahdc_track_phi->Write(h1_ahdc_track_phi->GetName());
        h1_ahdc_track_vz->Write(h1_ahdc_track_vz->GetName());

        h1_ahdc_track_delta_phi->Write(h1_ahdc_track_delta_phi->GetName());

        h1_elastics_ahdc_track_delta_phi->Write(h1_elastics_ahdc_track_delta_phi->GetName());
        h1_elastics_ahdc_track_nhits->Write(h1_elastics_ahdc_track_nhits->GetName());
        h1_elastics_ahdc_track_chi2->Write(h1_elastics_ahdc_track_chi2->GetName());
        h1_elastics_ahdc_track_p->Write(h1_elastics_ahdc_track_p->GetName());
        h1_elastics_ahdc_track_theta->Write(h1_elastics_ahdc_track_theta->GetName());
        h1_elastics_ahdc_track_phi->Write(h1_elastics_ahdc_track_phi->GetName());
        h1_elastics_ahdc_track_vz->Write(h1_elastics_ahdc_track_vz->GetName());
        h2_elastics_ahdc_track_p_vs_dEdx->Write(h2_elastics_ahdc_track_p_vs_dEdx->GetName());
        h1_elastics_electron_p->Write(h1_elastics_electron_p->GetName());
        h1_elastics_electron_theta->Write(h1_elastics_electron_theta->GetName());
        h1_elastics_electron_phi->Write(h1_elastics_electron_phi->GetName());
        h1_elastics_electron_vz->Write(h1_elastics_electron_vz->GetName());

        h1_elastics_delta_vz->Write(h1_elastics_delta_vz->GetName());
        h1_elastics_expected_track_p->Write(h1_elastics_expected_track_p->GetName());
        h1_elastics_expected_track_theta->Write(h1_elastics_expected_track_theta->GetName());
        h1_elastics_expected_track_phi->Write(h1_elastics_expected_track_phi->GetName());

        h1_elastics_diff_p->Write(h1_elastics_diff_p->GetName());
        h1_elastics_diff_theta->Write(h1_elastics_diff_theta->GetName());
        h1_elastics_diff_phi->Write(h1_elastics_diff_phi->GetName());

    }


};


/**
 * @brief Progress bar
 * 
 * @param state is a number below bar_length
 * @param bar_length is 100 per default
 */
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


#endif
