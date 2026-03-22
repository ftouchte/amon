#ifndef KF_MON_DATA_H
#define KF_MON_DATA_H

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"

struct Particle_kinematics {
    int row = -1;
    double px = -999, py = -999, pz = -999;
    double vz = -999;
    double p = -999, theta = -999, phi = -999;
};

// To be added in a library
void progressBar(int state, int bar_length = 100);
// Dispaly h_sel above h_all, useful to see before and after cuts
TCanvas* superimpose_histograms(TH1D* h_all, TH1D* h_sel, const char * name);
// pT versus Sum_ADC, dispaly the elastics cuts, just to save space
TCanvas* display_elastics_cuts(TH2D* h, const char * name);
// Try to fit histograms with 2 gaussians
TCanvas* fit_histogram(TH1D* h, const char * name);
// Small study to count the number of (AHDC, ATOF) matching
void count_atof_matching(std::string filename, TFile * rootFile);
// special routine for simulation
void run_simulation(std::string filename, TFile * rootFile);
// extract error from 2D histograms : the flag is to distinguish residual versus ADC and Time
//std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name);
//std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name, bool flag = false);
std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name, bool flag = false, TDirectory * dir = nullptr);
bool compareWedge(int s0, int l0, int w0, int sector, int layer, int wedge);
bool IsAMatch(int s0, int l0, int w0, int sector, int layer, int wedge);
int layer2number(int digit);
int slc2wire(int sector, int layer, int component);
void wire2slc(int wire, int & sector, int & layer, int & component);

template <typename T>
void writeHistoVector(std::vector<T*> & histos, TDirectory* dir, const char * title) {
    dir->cd();
    TCanvas* c1 = new TCanvas(TString::Format("c_%s", title).Data(),"");
    c1->Divide(3,3);
    for (int i = 0; i < (int) histos.size(); i++) {
        c1->cd(i+1);
        histos[i]->Draw();
    }
    c1->Write(title);
}

/*
class Histograms {

public:
    Histograms();
    ~Histograms();

public:
    
    TH1D* H1_cuts;
    std::vector<std::string> cutNames = {"trigger electron", "w2cut & 1st electron", "all tracks", "tracks n_hits >= 6", "best track (dphi)", "delta phi cut", "deuteron", "proton"};
    // electrons
    TH1D* H1_W2_all;
    TH1D* H1_W2_sel;
    TH1D* H1_Q2_all;
    TH1D* H1_Q2_sel;
    TH1D* H1_nu_all;
    TH1D* H1_nu_sel;
    TH1D* H1_vz;
    // ahdc track
    TH1D* H1_delta_phi_nosym; 
    TH1D* H1_delta_phi_nosym_sel; 
    TH1D* H1_delta_phi_sym; 
    TH1D* H1_delta_phi_sym_sel;
    TH2D* H2_corr_pTe_Sadc;
    TH2D* H2_corr_p_Sadc;
    TH2D* H2_corr_p_dEdx;
    TH2D* H2_corr_time_adc;
    // correlations electron vs ahdc_track
    std::vector<TH2D*> H2_corr_phi;
    std::vector<TH2D*> H2_corr_pT;
    std::vector<TH2D*> H2_corr_vz;
    std::vector<TH1D*> H1_delta_vz;
    std::vector<TH1D*> H1_delta_phi;
    // kinematics 1D
    std::vector<TH1D*> H1_track_pT;
    std::vector<TH1D*> H1_track_theta;
    std::vector<TH1D*> H1_track_phi;
    std::vector<TH1D*> H1_track_residual;
    std::vector<TH1D*> H1_track_residual_LR;
    std::vector<TH1D*> H1_track_chi2;
    std::vector<TH1D*> H1_track_nhits;
    std::vector<TH1D*> H1_track_sum_residual;
    std::vector<TH1D*> H1_electron_pT;
    std::vector<TH1D*> H1_electron_theta;
    std::vector<TH1D*> H1_electron_phi;
    // for elastics only, comparison to expected track kinematics given the scattered electron kinematics
    std::vector<TH1D*> H1_diff_pT;
    std::vector<TH1D*> H1_diff_theta;
    std::vector<TH1D*> H1_diff_phi;
    std::vector<TH2D*> H2_corr_residual_time;
    std::vector<TH2D*> H2_corr_residual_ADC;
    std::vector<TH2D*> H2_corr_residual_vz;
    std::vector<TH1D*> H1_residual_per_slayer;
    std::vector<TH2D*> H2_corr_residual_per_slayer_vz;
    std::vector<TH2D*> H2_time2distance;
    std::vector<TH1D*> H1_time;
    std::vector<TH1D*> H1_timeOverThreshold;
    std::vector<TH1D*> H1_amplitude;
    std::vector<TH1D*> H1_distance;
    // ATOF
    TH2D* H2_corr_atof_wedge_vz;
    TH1D* H1_diff_atof_wedge_vz;
    TH1D* H1_nmatched_bar;
    TH2D* H2_corr_atof_bar_vz;
    TH1D* H1_diff_atof_bar_vz;
    TH1D* H1_residual_with_atof;
    TH1D* H1_residual_without_atof;
    TH1D* H1_delta_phi_with_atof;
    TH1D* H1_delta_phi_without_atof;
    TH1D* H1_delta_vz_with_atof;
    TH1D* H1_delta_vz_without_atof;
    TH1D* H1_track_theta_with_atof;
    TH1D* H1_track_theta_without_atof;
    TH1D* H1_track_atof_region;

    TH1D* H1_track_atof_s2_sigma_z;
    TH1D* H1_track_atof_s2_sigma_phi;

};


Histograms::Histograms() {
    H1_cuts = new TH1D("cuts", "Nb. events over cuts", 8, 0, 8);
    // std::vector<std::string> cutNames = {"trigger electron", "w2cut & 1st electron", "all tracks", "tracks n_hits >= 6", "best track (dphi)", "delta phi cut", "deuteron", "proton"};
    // electrons
    H1_W2_all = new TH1D("W2", "W^{2}; W^{2} (GeV^{2}); count", 100, 3.2, 7);
    H1_W2_sel = new TH1D("W2_sel", "W^{2}; W^{2} (GeV^{2}); count", 100, 3.2, 7);
    H1_Q2_all = new TH1D("Q2", "Q^{2}; Q^{2} (GeV^{2}); count", 50, 0, 0.4);
    H1_Q2_sel = new TH1D("Q2_sel", "Q^{2}; Q^{2} (GeV^{2}); count", 50, 0, 0.4);
    H1_nu_all = new TH1D("nu", "#nu = E - E'; #nu (GeV); count", 50, 0, 2.24);
    H1_nu_sel = new TH1D("nu_sel", "#nu = E - E'; #nu (GeV); count", 50, 0, 2.24);
    H1_vz = new TH1D("vz_electron", "vz; vz (cm); count", 50, -60, 40);
    // ahdc track
    H1_delta_phi_nosym = new TH1D("delta_phi_nosym", "#Delta #phi; #Delta #phi (deg); count", 100, -360, 360); 
    H1_delta_phi_nosym_sel = new TH1D("delta_phi_nosym_sel", "#Delta #phi; #Delta #phi (deg); count", 100, -360, 360); 
    H1_delta_phi_sym = new TH1D("delta_phi_sym", "#Delta #phi; #Delta #phi (deg); count", 100, -90, 90); 
    H1_delta_phi_sym_sel = new TH1D("delta_phi_sym_sel", "#Delta #phi; #Delta #phi (deg); count", 100, -90, 90);
    H2_corr_pTe_Sadc = new TH2D("corr_pTe_Sadc", "#Sigma ADC vs pT; pT (MeV); #Sigma ADC", 50, 190, 400, 50, 0, 15000);
    H2_corr_p_Sadc = new TH2D("corr_p_Sadc", "#Sigma ADC vs p; p (MeV); #Sigma ADC", 50, 80, 600, 50, 0, 15000);
    H2_corr_p_dEdx = new TH2D("corr_p_dEdx", "dEdx vs p; p (MeV); dEdx (MeV/mm)", 50, 80, 600, 50, 0, 180);
    H2_corr_time_adc = new TH2D("corr_time_adc", "time vs ADC; ADC; time", 50, 0, 3700, 50, 0, 250);
    // correlations electron vs ahdc_track
    // std::vector<TH2D*> H2_corr_phi;
        H2_corr_phi.push_back(new TH2D("corr_phi_all", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 100, 0, 361, 50, 0, 361)); // all elastics
        H2_corr_phi.push_back(new TH2D("corr_phi_deuteron", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 100, 0, 361, 50, 0, 361)); // deuteron
        H2_corr_phi.push_back(new TH2D("corr_phi_proton", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 100, 0, 361, 50, 0, 361)); // proton
    // std::vector<TH2D*> H2_corr_pT;
        H2_corr_pT.push_back(new TH2D("corr_pT_all", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 50, 80, 400)); // all elastics
        H2_corr_pT.push_back(new TH2D("corr_pT_deuteron", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 50, 80, 400)); // deuteron
        H2_corr_pT.push_back(new TH2D("corr_pT_proton", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 50, 80, 400)); // proton
    // std::vector<TH2D*> H2_corr_vz;
        H2_corr_vz.push_back(new TH2D("corr_vz_all", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 50, -20, 16)); // all elastics
        H2_corr_vz.push_back(new TH2D("corr_vz_deuteron", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 50, -20, 16)); // deuteron
        H2_corr_vz.push_back(new TH2D("corr_vz_proton", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 50, -20, 16)); // proton
    // std::vector<TH1D*> H1_delta_vz;
        H1_delta_vz.push_back(new TH1D("delta_vz_all", "#Delta vz = vz^{(electron)} - vz^{(AHDC track)}; #Delta vz (cm); #count", 100, -20, 10)); // all elastics
        H1_delta_vz.push_back(new TH1D("delta_vz_deuteron", "#Delta vz = vz^{(electron)} - vz^{(AHDC track)}; #Delta vz (cm); #count", 100, -20, 10)); // deuteron
        H1_delta_vz.push_back(new TH1D("delta_vz_proton", "#Delta vz = vz^{(electron)} - vz^{(AHDC track)}; #Delta vz (cm); #count", 100, -20, 10)); // proton
    // std::vector<TH1D*> H1_delta_phi;
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_all", "#Delta #phi; #Delta #phi (deg); count", 100, -1.05*delta_phi_width, 1.05*delta_phi_width)); // all elastics
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_deuteron", "#Delta #phi; #Delta #phi (deg); count", 100, -1.05*delta_phi_width, 1.05*delta_phi_width)); // deuteron
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_proton", "#Delta #phi; #Delta #phi (deg); count", 100, -1.05*delta_phi_width, 1.05*delta_phi_width)); // proton
    // kinematics 1D
    // std::vector<TH1D*> H1_track_pT;
        H1_track_pT.push_back(new TH1D("track_pT_all", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // all elastics
        H1_track_pT.push_back(new TH1D("track_pT_deuteron", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // deuteron
        H1_track_pT.push_back(new TH1D("track_pT_proton", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // proton
    // std::vector<TH1D*> H1_track_theta;
        H1_track_theta.push_back(new TH1D("track_theta_all", "#theta ; #theta (deg); count", 50, 0, 180)); // all elastics
        H1_track_theta.push_back(new TH1D("track_theta_deuteron", "#theta ; #theta (deg); count", 50, 0, 180)); // deuteron
        H1_track_theta.push_back(new TH1D("track_theta_proton", "#theta ; #theta (deg); count", 50, 0, 180)); // proton
    // std::vector<TH1D*> H1_track_phi;
        H1_track_phi.push_back(new TH1D("track_phi_all", "#phi_{t} ; #phi_{t} (deg); count", 50, 0, 361)); // all elastics
        H1_track_phi.push_back(new TH1D("track_phi_deuteron", "#phi_{t} ; #phi_{t} (deg); count", 50, 0, 361)); // deuteron
        H1_track_phi.push_back(new TH1D("track_phi_proton", "#phi_{t} ; #phi_{t} (deg); count", 50, 0, 361)); // proton
    // std::vector<TH1D*> H1_track_residual;
        H1_track_residual.push_back(new TH1D("track_residual_all", "residual ; residual (mm); count", 100, -3, 3)); // all elastics
        H1_track_residual.push_back(new TH1D("track_residual_deuteron", "residual ; residual (mm); count", 100, -3, 3)); // deuteron
        H1_track_residual.push_back(new TH1D("track_residual_proton", "residual ; residual (mm); count", 100, -3, 3)); // proton
    // std::vector<TH1D*> H1_track_residual_LR;
        H1_track_residual_LR.push_back(new TH1D("track_residual_LR_all", "residual_LR ; residual_LR (mm); count", 50, -3, 3)); // all elastics
        H1_track_residual_LR.push_back(new TH1D("track_residual_LR_deuteron", "residual_LR ; residual_LR (mm); count", 50, -3, 3)); // deuteron
        H1_track_residual_LR.push_back(new TH1D("track_residual_LR_proton", "residual_LR ; residual_LR (mm); count", 50, -3, 3)); // proton
    // std::vector<TH1D*> H1_track_chi2;
        H1_track_chi2.push_back(new TH1D("track_chi2_all", "chi2 ; chi2; count", 50, 0, 5)); // all elastics
        H1_track_chi2.push_back(new TH1D("track_chi2_deuteron", "chi2 ; chi2; count", 50, 0, 5)); // deuteron
        H1_track_chi2.push_back(new TH1D("track_chi2_proton", "chi2 ; chi2; count", 50, 0, 5)); // proton
    // std::vector<TH1D*> H1_track_nhits;
        H1_track_nhits.push_back(new TH1D("track_nhits_all", "number of hits per track; nhits; count", 12, 5, 17)); // all elastics
        H1_track_nhits.push_back(new TH1D("track_nhits_deuteron", "number of hits per track; nhits; count", 12, 5, 17)); // deuteron
        H1_track_nhits.push_back(new TH1D("track_nhits_proton", "number of hits per track; nhits; count", 12, 5, 17)); // proton
    // std::vector<TH1D*> H1_track_sum_residual;
        H1_track_sum_residual.push_back(new TH1D("track_sum_residual_all", "Sum of residuals per track; sum_residual (mm); count", 50, -10, 5)); // all elastics
        H1_track_sum_residual.push_back(new TH1D("track_sum_residual_deuteron", "Sum of residuals per track; sum_residual (mm); count", 50, -10, 5)); // deuteron
        H1_track_sum_residual.push_back(new TH1D("track_sum_residual_proton", "Sum of residuals per track; sum_residual (mm); count", 50, -10, 5)); // proton
    // std::vector<TH1D*> H1_electron_pT;
        H1_electron_pT.push_back(new TH1D("electron_pT_all", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // all elastics
        H1_electron_pT.push_back(new TH1D("electron_pT_deuteron", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // deuteron
        H1_electron_pT.push_back(new TH1D("electron_pT_proton", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // proton
    // std::vector<TH1D*> H1_electron_theta;
        H1_electron_theta.push_back(new TH1D("electron_theta_all", "#theta_{e} ; #theta_{e} (deg); count", 50, 0, 20)); // all elastics
        H1_electron_theta.push_back(new TH1D("electron_theta_deuteron", "#theta_{e} ; #theta_{e} (deg); count", 50, 0, 20)); // deuteron 
        H1_electron_theta.push_back(new TH1D("electron_theta_proton", "#theta_{e} ; #theta_{e} (deg); count", 50, 0, 20)); // proton
    // std::vector<TH1D*> H1_electron_phi;
        H1_electron_phi.push_back(new TH1D("electron_phi_all", "#phi_{e} ; #phi_{e} (deg); count", 50, 0, 361)); // all elastics
        H1_electron_phi.push_back(new TH1D("electron_phi_deuteron", "#phi_{e} ; #phi_{e} (deg); count", 50, 0, 361)); // deuteron
        H1_electron_phi.push_back(new TH1D("electron_phi_proton", "#phi_{e} ; #phi_{e} (deg); count", 50, 0, 361)); // proton
    // for elastics only, comparison to expected track kinematics given the scattered electron kinematics
    // std::vector<TH1D*> H1_diff_pT;
        H1_diff_pT.push_back(new TH1D("diff_pT_all", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 400));
        H1_diff_pT.push_back(new TH1D("diff_pT_deuteron", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 140));
        H1_diff_pT.push_back(new TH1D("diff_pT_proton", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 400));
    // std::vector<TH1D*> H1_diff_theta;
        H1_diff_theta.push_back(new TH1D("diff_theta_all", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -80, 80));
        H1_diff_theta.push_back(new TH1D("diff_theta_deuteron", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -80, 80));
        H1_diff_theta.push_back(new TH1D("diff_theta_proton", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -80, 80));
    // std::vector<TH1D*> H1_diff_phi;
        H1_diff_phi.push_back(new TH1D("diff_phi_all", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));
        H1_diff_phi.push_back(new TH1D("diff_phi_deuteron", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));
        H1_diff_phi.push_back(new TH1D("diff_phi_proton", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));
    // std::vector<TH2D*> H2_corr_residual_time;
        H2_corr_residual_time.push_back(new TH2D("corr_residual_time_all", "residual vs time; time (ns); residual (mm)", 50, 0, 250, 50, -3, 3)); // all elastics
        H2_corr_residual_time.push_back(new TH2D("corr_residual_time_deuteron", "residual vs time; time (ns); residual (mm)", 50, 0, 250, 50, -3, 3)); // deuteron
        H2_corr_residual_time.push_back(new TH2D("corr_residual_time_proton", "residual vs time; time (ns); residual (mm)", 50, 0, 250, 50, -3, 3)); // proton
    // std::vector<TH2D*> H2_corr_residual_ADC;
        H2_corr_residual_ADC.push_back(new TH2D("corr_residual_ADC_all", "residual vs ADC; ADC; residual (mm)", 50, 0, 3700, 50, -3, 3)); // all elastics
        H2_corr_residual_ADC.push_back(new TH2D("corr_residual_ADC_deuteron", "residual vs ADC; ADC; residual (mm)", 50, 0, 3700, 50, -3, 3)); // deuteron
        H2_corr_residual_ADC.push_back(new TH2D("corr_residual_ADC_proton", "residual vs ADC; ADC; residual (mm)", 50, 0, 3700, 50, -3, 3)); // proton
    // std::vector<TH2D*> H2_corr_residual_vz;
        H2_corr_residual_vz.push_back(new TH2D("corr_residual_vz_all", "residual vs vz; vz (cm); residual (mm)", 50, -25, 20, 50, -3, 3)); // all elastics
        H2_corr_residual_vz.push_back(new TH2D("corr_residual_vz_deuteron", "residual vs vz (cm); vz; residual (mm)", 50, -25, 20, 50, -3, 3)); // deuteron
        H2_corr_residual_vz.push_back(new TH2D("corr_residual_vz_proton", "residual vs vz; vz (cm); residual (mm)", 50, -25, 30, 20, -3, 3)); // proton
    // std::vector<TH1D*> H1_residual_per_slayer;
        H1_residual_per_slayer.push_back(new TH1D("track_residual_slayer_1", "residual per super layer 1; residual (mm); count", 50, -3, 3));
        H1_residual_per_slayer.push_back(new TH1D("track_residual_slayer_2", "residual per super layer 2; residual (mm); count", 50, -3, 3));
        H1_residual_per_slayer.push_back(new TH1D("track_residual_slayer_3", "residual per super layer 3; residual (mm); count", 50, -3, 3));
        H1_residual_per_slayer.push_back(new TH1D("track_residual_slayer_4", "residual per super layer 4; residual (mm); count", 50, -3, 3));
        H1_residual_per_slayer.push_back(new TH1D("track_residual_slayer_5", "residual per super layer 5; residual (mm); count", 50, -3, 3));
    // std::vector<TH2D*> H2_corr_residual_per_slayer_vz;
        H2_corr_residual_per_slayer_vz.push_back(new TH2D("track_residual_slayer_vz_1", "residual per super layer 1; vz (cm); residual (mm)",  50, -25, 20, 50, -3, 3));
        H2_corr_residual_per_slayer_vz.push_back(new TH2D("track_residual_slayer_vz_2", "residual per super layer 2; vz (cm); residual (mm)",  50, -25, 20, 50, -3, 3));
        H2_corr_residual_per_slayer_vz.push_back(new TH2D("track_residual_slayer_vz_3", "residual per super layer 3; vz (cm); residual (mm)",  50, -25, 20, 50, -3, 3));
        H2_corr_residual_per_slayer_vz.push_back(new TH2D("track_residual_slayer_vz_4", "residual per super layer 4; vz (cm); residual (mm)",  50, -25, 20, 50, -3, 3));
        H2_corr_residual_per_slayer_vz.push_back(new TH2D("track_residual_slayer_vz_5", "residual per super layer 5; vz (cm); residual (mm)",  50, -25, 20, 50, -3, 3));
    // std::vector<TH2D*> H2_time2distance;
        H2_time2distance.push_back(new TH2D("corr_t2d_all", "time2distance;time (ns); distance (mm)", 100, 0, 320, 100, 0, 4));
        H2_time2distance.push_back(new TH2D("corr_t2d_deuteron", "time2distance;time (ns); distance (mm)", 100, 0, 320, 100, 0, 4));
        H2_time2distance.push_back(new TH2D("corr_t2d_proton", "time2distance;time (ns); distance (mm)", 100, 0, 320, 100, 0, 4));
    // std::vector<TH1D*> H1_time;
        H1_time.push_back(new TH1D("time_all", "time ; time (ns); count",  100, 0, 320)); // all elastics
        H1_time.push_back(new TH1D("time_deuteron", "time ; time (ns); count",  100, 0, 320)); // deuteron 
        H1_time.push_back(new TH1D("time_proton", "time ; time (ns); count",  100, 0, 320)); // proton
    // std::vector<TH1D*> H1_timeOverThreshold;
        H1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_all", "timeOverThreshold ; timeOverThreshold (ns); count",  100, 250, 700)); // all elastics
        H1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_deuteron", "timeOverThreshold ; timeOverThreshold (ns); count",  100, 250, 700)); // deuteron 
        H1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_proton", "timeOverThreshold ; timeOverThreshold (ns); count",  100, 250, 700)); // proton
    // std::vector<TH1D*> H1_amplitude;
        H1_amplitude.push_back(new TH1D("amplitude_all", "amplitude ; amplitude (ADC); count",  100, 0, 4000)); // all elastics
        H1_amplitude.push_back(new TH1D("amplitude_deuteron", "amplitude ; amplitude (ADC); count",  100, 0, 4000)); // deuteron 
        H1_amplitude.push_back(new TH1D("amplitude_proton", "amplitude ; amplitude (ADC); count",  100, 0, 4000)); // proton
    // std::vector<TH1D*> H1_distance;
        H1_distance.push_back(new TH1D("AHDC::hits:doca_all", "AHDC::hits:doca ; AHDC::hits:doca (mm); count",  50, 0, 4)); // all elastics
        H1_distance.push_back(new TH1D("AHDC::hits:doca_deuteron", "AHDC::hits:doca ; AHDC::hits:doca (mm); count",  50, 0, 4)); // deuteron 
        H1_distance.push_back(new TH1D("AHDC::hits:doca_proton", "AHDC::hits:doca ; AHDC::hits:doca (mm); count",  50, 0, 4)); // proton
    // ATOF
    H2_corr_atof_wedge_vz = new TH2D("corr_atof_wedge_vz", "ATOF wedge vz versus electron vz;electron vz (cm); ATOF wedge vz (cm)", 50, -30, 20, 10, -14, 14);
    H1_diff_atof_wedge_vz = new TH1D("diff_atof_wedge_vz", "#Delta vz = vz^{(electron)} - vz^{(ATOF wedge)}; #Delta vz (cm); count", 50, -30, 10);
    H1_nmatched_bar = new TH1D("nmatched_bar", "number of bars for a matched wegde; nbars; count", 6, 0, 6);
    H2_corr_atof_bar_vz = new TH2D("corr_atof_vz_all", "ATOF bar vz versus electron vz;electron vz (cm); ATOF bar vz (cm)", 50, -30, 20, 100, -20, 20);
    H1_diff_atof_bar_vz = new TH1D("diff_atof_vz_all", "#Delta vz = vz^{(electron)} - vz^{(ATOF bar)}; #Delta vz (cm); count", 50, -30, 10);
    H1_residual_with_atof = new TH1D("residual_with_atof", "residual; residual (mm); count", 50, -3, 3);
    H1_residual_without_atof = new TH1D("residual_without_atof", "residual; residual (mm); count", 50, -3, 3);
    H1_delta_phi_with_atof = new TH1D("delta_phi_with_atof", "#Delta #phi; #Delta #phi (deg); count", 100, -1.05*delta_phi_width, 1.05*delta_phi_width);
    H1_delta_phi_without_atof = new TH1D("delta_phi_without_atof", "#Delta #phi; #Delta #phi (deg); count", 100, -1.05*delta_phi_width, 1.05*delta_phi_width);
    H1_delta_vz_with_atof = new TH1D("delta_vz_with_atof", "#Delta vz = vz^{(electron)} - vz^{(AHDC track)}; #Delta vz (cm); #count", 100, -20, 10);
    H1_delta_vz_without_atof = new TH1D("delta_vz_without_atof", "#Delta vz = vz^{(electron)} - vz^{(AHDC track)}; #Delta vz (cm); #count", 100, -20, 10);
    H1_track_theta_with_atof = new TH1D("track_theta_with_atof", "#theta ; #theta (deg); count", 50, 0, 180);
    H1_track_theta_without_atof = new TH1D("track_theta_without_atof", "#theta ; #theta (deg); count", 50, 0, 180);
    H1_track_atof_region = new TH1D("track_atof_region", "ATOF region ; region; count", 4, 0, 4);

    H1_track_atof_s2_sigma_z = new TH1D("H1_track_atof_s2_sigma_z", "#sigma_{z} ; #sigma_{z} (cm); count", 100, 2.95, 3.8);
    H1_track_atof_s2_sigma_phi = new TH1D("H1_track_atof_s2_sigma_phi", "#sigma_{#phi} ; #sigma_{#phi} (deg); count", 100, 4.9, 8);
}


Histograms::~Histograms() {

}*/








#endif
