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
// Given a list of rows, fill the corresponding histograms and entries reading the bank
void fill_histogram(std::vector<TH1D*> histos, std::vector<const char *> entries, hipo::bank & bank, std::vector<int> & rows);
// Small study to count the number of (AHDC, ATOF) matching
void count_atof_matching(std::string filename, TFile * rootFile);
// special routine for simulation
void run_simulation(std::string filename, TFile * rootFile);
// extract error from 2D histograms : the flag is to distinguish residual versus ADC and Time
//std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name);
//std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name, bool flag = false);
std::pair<TCanvas*, TGraph*> projectionY(TH2D* h, int nbins_per_projection, const char * name, bool flag = false, TDirectory * dir = nullptr);

#endif