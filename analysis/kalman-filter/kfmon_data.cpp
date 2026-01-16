/***********************************************
 * Monitoring du Kalman Filter (case of real 
 * data)
 *
 * @author Felix Touchte Codjo
 * @date January 07, 2026
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
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "Math/PdfFuncMathCore.h"
#include "TText.h"
#include "THStack.h"
#include "TArrow.h"

#include "futils.h"

// utilities
void progressBar(int state, int bar_length = 100);
//////////////////////////////

// define some constants
const double beam_energy = 2.23951; // incident energy if the electron, GeV
const double electron_mass = 0.511e-3; // energy mass of electron, GeV
const double proton_mass = 938.272e-3; // energy mass of proton, GeV
const double helium_mass = 3.73; // energy mass of Helium-4, GeV
const double deuteron_mass = 1.875; // energy mass of Deuterium, GeV

const double delta_phi_width = 45;

struct Particle_kinematics {
    int row = -1;
    double px = -999, py = -999, pz = -999;
    double vz = -999;
    double p = -999, theta = -999, phi = -999;
};

TCanvas* superimpose_histograms(TH1D* h_all, TH1D* h_sel, const char * name);
TCanvas* display_elastics_cuts(TH2D* h, const char * name);
TCanvas* fit_histogram(TH1D* h, const char * name);
void fill_histogram(std::vector<TH1D*> histos, std::vector<const char *> entries, hipo::bank & bank, std::vector<int> & rows);
void count_atof_matching(std::string filename, TFile * rootFile);

// main routine
int main(int argc, char const *argv[]) {
    // record start time
    auto start = std::chrono::high_resolution_clock::now();

    std::string filename, output;

    if (argc < 2) {
        printf("Please provide 1 file name and 1 output name...\n");
        // return error code
        return 1;
    }
    
    // set filename and initialise reader
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-22712.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-realrun-22712-p0v6.hipo";

    filename = argv[1];
    output = argv[2];

    //const char * output = "./output/kfmon_data.root";
    printf("> filename : %s\n", filename.c_str());
    hipo::reader  reader(filename.c_str());
    hipo::dictionary factory;
    reader.readDictionary(factory);

    // bank definition
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  recBank(factory.getSchema("REC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;

    // example of 1D histograms
    //TH1D* H1_mctime = new TH1D("mctime", "mctime; mctime (ns); count", 100, 0, 250); 
    // electrons
    TH1D* H1_W2_all = new TH1D("W2", "W^{2}; W^{2} (GeV^{2}); count", 50, 2, 7);
    TH1D* H1_W2_sel = new TH1D("W2_sel", "W^{2}; W^{2} (GeV^{2}); count", 50, 2, 7);
    TH1D* H1_Q2_all = new TH1D("Q2", "Q^{2}; Q^{2} (GeV^{2}); count", 50, 0, 0.4);
    TH1D* H1_Q2_sel = new TH1D("Q2_sel", "Q^{2}; Q^{2} (GeV^{2}); count", 50, 0, 0.4);
    TH1D* H1_nu_all = new TH1D("nu", "#nu = E - E'; #nu (GeV); count", 50, 0, 2.24);
    TH1D* H1_nu_sel = new TH1D("nu_sel", "#nu = E - E'; #nu (GeV); count", 50, 0, 2.24);
    TH1D* H1_vz = new TH1D("vz_electron", "vz; vz (cm); count", 50, -60, 40);
    // ahdc track
    TH1D* H1_delta_phi_nosym = new TH1D("delta_phi_nosym", "#Delta #phi; #Delta #phi (deg); count", 50, -360, 360); 
    TH1D* H1_delta_phi_nosym_sel = new TH1D("delta_phi_nosym_sel", "#Delta #phi; #Delta #phi (deg); count", 50, -360, 360); 
    TH1D* H1_delta_phi_sym = new TH1D("delta_phi_sym", "#Delta #phi; #Delta #phi (deg); count", 50, -180, 180); 
    TH1D* H1_delta_phi_sym_sel = new TH1D("delta_phi_sym_sel", "#Delta #phi; #Delta #phi (deg); count", 50, -180, 180);
    TH2D* H2_corr_pTe_Sadc = new TH2D("corr_pTe_Sadc", "#Sigma ADC vs pT_{e}; pT (MeV); #Sigma ADC", 50, 190, 400, 100, 0, 15000);
    TH2D* H2_corr_pTt_Sadc = new TH2D("corr_pTt_Sadc", "#Sigma ADC vs pT_{t}; pT (MeV); #Sigma ADC", 50, 80, 600, 100, 0, 15000);
    // correlations electron vs ahdc_track
    std::vector<TH2D*> H2_corr_phi;
        H2_corr_phi.push_back(new TH2D("corr_phi_all", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 50, 0, 361, 100, 0, 361)); // all elastics
        H2_corr_phi.push_back(new TH2D("corr_phi_deuteron", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 50, 0, 361, 100, 0, 361)); // deuteron
        H2_corr_phi.push_back(new TH2D("corr_phi_proton", "#phi_{t} vs #phi_{e}; #phi_{e} (deg); #phi_{t} (deg)", 50, 0, 361, 100, 0, 361)); // proton
    std::vector<TH2D*> H2_corr_pT;
        H2_corr_pT.push_back(new TH2D("corr_pT_all", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 100, 80, 400)); // all elastics
        H2_corr_pT.push_back(new TH2D("corr_pT_deuteron", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 100, 80, 400)); // deuteron
        H2_corr_pT.push_back(new TH2D("corr_pT_proton", "pT_{t} vs pT_{e}; pT_{e} (MeV); pT_{t} (MeV)", 50, 190, 400, 100, 80, 400)); // proton
    std::vector<TH2D*> H2_corr_vz;
        H2_corr_vz.push_back(new TH2D("corr_vz_all", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 100, -20, 16)); // all elastics
        H2_corr_vz.push_back(new TH2D("corr_vz_deuteron", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 100, -20, 16)); // deuteron
        H2_corr_vz.push_back(new TH2D("corr_vz_proton", "vz_{t} vs vz_{e}; vz_{e} (cm); vz_{t} (cm)", 50, -35, 16, 100, -20, 16)); // proton
    std::vector<TH1D*> H1_delta_vz;
        H1_delta_vz.push_back(new TH1D("delta_vz_all", "#Delta vz = vz_{e} - vz_{t}; #Delta vz (cm); #count", 50, -30, 30)); // all elastics
        H1_delta_vz.push_back(new TH1D("delta_vz_deuteron", "#Delta vz = vz_{e} - vz_{t}; #Delta vz (cm); #count", 50, -30, 30)); // deuteron
        H1_delta_vz.push_back(new TH1D("delta_vz_proton", "#Delta vz = vz_{e} - vz_{t}; #Delta vz (cm); #count", 50, -30, 30)); // proton
    std::vector<TH1D*> H1_delta_phi;
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_all", "#Delta #phi; #Delta #phi (deg); count", 50, -1.05*delta_phi_width, 1.05*delta_phi_width)); // all elastics
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_deuteron", "#Delta #phi; #Delta #phi (deg); count", 50, -1.05*delta_phi_width, 1.05*delta_phi_width)); // deuteron
        H1_delta_phi.push_back(new TH1D("delta_phi_sym_proton", "#Delta #phi; #Delta #phi (deg); count", 50, -1.05*delta_phi_width, 1.05*delta_phi_width)); // proton
    // kinematics 1D
    std::vector<TH1D*> H1_track_pT;
        H1_track_pT.push_back(new TH1D("track_pT_all", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // all elastics
        H1_track_pT.push_back(new TH1D("track_pT_deuteron", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // deuteron
        H1_track_pT.push_back(new TH1D("track_pT_proton", "pT_{t} ; pT_{t} (MeV); count", 50, 0, 1000)); // proton
    std::vector<TH1D*> H1_track_theta;
        H1_track_theta.push_back(new TH1D("track_theta_all", "#theta_{t} ; #theta_{t} (MeV); count", 50, 0, 180)); // all elastics
        H1_track_theta.push_back(new TH1D("track_theta_deuteron", "#theta_{t} ; #theta_{t} (MeV); count", 50, 0, 180)); // deuteron
        H1_track_theta.push_back(new TH1D("track_theta_proton", "#theta_{t} ; #theta_{t} (MeV); count", 50, 0, 180)); // proton
    std::vector<TH1D*> H1_track_phi;
        H1_track_phi.push_back(new TH1D("track_phi_all", "#phi_{t} ; #phi_{t} (MeV); count", 50, 0, 361)); // all elastics
        H1_track_phi.push_back(new TH1D("track_phi_deuteron", "#phi_{t} ; #phi_{t} (MeV); count", 50, 0, 361)); // deuteron
        H1_track_phi.push_back(new TH1D("track_phi_proton", "#phi_{t} ; #phi_{t} (MeV); count", 50, 0, 361)); // proton
    std::vector<TH1D*> H1_track_residual;
        H1_track_residual.push_back(new TH1D("track_residual_all", "residual ; residual (mm); count", 50, -3, 3)); // all elastics
        H1_track_residual.push_back(new TH1D("track_residual_deuteron", "residual ; residual (mm); count", 50, -3, 3)); // deuteron
        H1_track_residual.push_back(new TH1D("track_residual_proton", "residual ; residual (mm); count", 50, -3, 3)); // proton
    std::vector<TH1D*> H1_electron_pT;
        H1_electron_pT.push_back(new TH1D("electron_pT_all", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // all elastics
        H1_electron_pT.push_back(new TH1D("electron_pT_deuteron", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // deuteron
        H1_electron_pT.push_back(new TH1D("electron_pT_proton", "pT_{e} ; pT_{e} (MeV); count", 50, 190, 400)); // proton
    std::vector<TH1D*> H1_electron_theta;
        H1_electron_theta.push_back(new TH1D("electron_theta_all", "#theta_{e} ; #theta_{e} (MeV); count", 50, 0, 20)); // all elastics
        H1_electron_theta.push_back(new TH1D("electron_theta_deuteron", "#theta_{e} ; #theta_{e} (MeV); count", 50, 0, 20)); // deuteron 
        H1_electron_theta.push_back(new TH1D("electron_theta_proton", "#theta_{e} ; #theta_{e} (MeV); count", 50, 0, 20)); // proton
    std::vector<TH1D*> H1_electron_phi;
        H1_electron_phi.push_back(new TH1D("electron_phi_all", "#phi_{e} ; #phi_{e} (MeV); count", 50, 0, 361)); // all elastics
        H1_electron_phi.push_back(new TH1D("electron_phi_deuteron", "#phi_{e} ; #phi_{e} (MeV); count", 50, 0, 361)); // deuteron
        H1_electron_phi.push_back(new TH1D("electron_phi_proton", "#phi_{e} ; #phi_{e} (MeV); count", 50, 0, 361)); // proton
    // for elastics only, comparison to expected track kinematics given the scattered electron kinematics
    std::vector<TH1D*> H1_diff_pT;
        H1_diff_pT.push_back(new TH1D("diff_pT_all", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 400));
        H1_diff_pT.push_back(new TH1D("diff_pT_deuteron", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 140));
        H1_diff_pT.push_back(new TH1D("diff_pT_proton", "(track) pT_{recon} - pT_{expected} ; pT_{recon} - pT_{expected} (MeV); count", 50, -400, 400));
    std::vector<TH1D*> H1_diff_theta;
        H1_diff_theta.push_back(new TH1D("diff_theta_all", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -180, 100));
        H1_diff_theta.push_back(new TH1D("diff_theta_deuteron", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -180, 100));
        H1_diff_theta.push_back(new TH1D("diff_theta_proton", "(track) theta_{recon} - theta_{expected} ; theta_{recon} - theta_{expected} (deg); count", 50, -180, 100));
    std::vector<TH1D*> H1_diff_phi;
        H1_diff_phi.push_back(new TH1D("diff_phi_all", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));
        H1_diff_phi.push_back(new TH1D("diff_phi_deuteron", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));
        H1_diff_phi.push_back(new TH1D("diff_phi_proton", "(track) phi_{recon} - phi_{expected} ; phi_{recon} - phi_{expected} (deg); count", 50, -50, 50));

    // example of 2D histograms
    
    /////////////////////////
    // Loop over events
    /////////////////////////
    while( reader.next()){
        nevents++;
        // display progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        // load bank content for this event
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(trackBank);
        event.getStructure(hitBank);
        event.getStructure(recBank);

        // find the elastic electron
        //int electron_row = -1;
        Particle_kinematics electron; 
        Particle_kinematics expected_track;
        for (int recRow = 0; recRow < recBank.getRows(); recRow++) {
            // electron pid
            if (recBank.getInt("pid", recRow) == 11 && recBank.getShort("status", recRow) < 0) {
                // compute kinematic variables
                double px = recBank.getDouble("px", recRow);
                double py = recBank.getDouble("py", recRow);
                double pz = recBank.getDouble("pz", recRow);
                double p, theta, phi;
                futils::cart2polarDEG(px,py,pz,p, theta, phi);
                double scattered_beam_energy = sqrt(pow(p,2) + pow(electron_mass,2));
                double nu = beam_energy - scattered_beam_energy;
                double Q2 = 4*beam_energy*scattered_beam_energy*pow(sin((theta/2)*M_PI/180),2);
                double W2 = pow(deuteron_mass,2) + 2*deuteron_mass*nu - Q2;
                //double xB = Q2/(2*proton_mass*nu);
                H1_nu_all->Fill(nu);
                H1_Q2_all->Fill(Q2);
                H1_W2_all->Fill(W2);
                H1_vz->Fill(recBank.getDouble("vz", recRow));
                
                // we only look at the first trigger electron (i.e status < 0)
                // + cut on W2 to select elastics
                if (electron.row < 0 && W2 > 3.4 && W2 < 3.9) {
                        electron.row = recRow;
                        electron.px = px; // GeV
                        electron.py = py; // GeV
                        electron.pz = pz; // GeV
                        electron.vz = recBank.getDouble("vz", recRow); // cm
                        electron.p = p;
                        electron.theta = theta;
                        electron.phi = phi;
                        H1_nu_sel->Fill(nu);
                        H1_Q2_sel->Fill(Q2);
                        H1_W2_sel->Fill(W2);
                        // for elastics only, comparison to expected track kinematics given the scattered electron kinematics
                        expected_track.px = -scattered_beam_energy*sin(electron.theta*M_PI/180)*cos(electron.phi*M_PI/180)*1000; // MeV
                        expected_track.py = -scattered_beam_energy*sin(electron.theta*M_PI/180)*sin(electron.phi*M_PI/180)*1000; // MeV
                        expected_track.pz = (beam_energy - scattered_beam_energy*cos(electron.theta*M_PI/180))*1000; // MeV
                        expected_track.vz = electron.vz;
                        expected_track.row = electron.row;
                        futils::cart2polarDEG(expected_track.px, expected_track.py, expected_track.pz, expected_track.p, expected_track.theta, expected_track.phi);
                }
            }
        }

        // we have an elastic electron
        // we now look at the AHDC::kftrack
        Particle_kinematics ahdc_track;
        double dphi = 1e10; // arbitrary big number
        if (electron.row >= 0) {
            for (int trackRow = 0; trackRow < trackBank.getRows(); trackRow++) {
                // read kinematic variables
                double px = trackBank.getDouble("px", trackRow);
                double py = trackBank.getDouble("py", trackRow);
                double pz = trackBank.getDouble("pz", trackRow);
                double p, theta, phi;
                futils::cart2polarDEG(px,py,pz,p, theta, phi);
                double delta_phi = fabs(electron.phi - phi) - 180;
                // select the best track (the one with the best delta_phi)
                if (fabs(delta_phi) < dphi) {
                    ahdc_track.row = trackRow;
                    ahdc_track.px = px; // MeV
                    ahdc_track.py = py; // MeV
                    ahdc_track.pz = pz; // MeV
                    ahdc_track.vz = trackBank.getFloat("z", trackRow); // mm
                    ahdc_track.p = p;
                    ahdc_track.theta = theta;
                    ahdc_track.phi = phi;
                    dphi = delta_phi;
                }
            }
        }
        // at this stage we have an electron and the best ahdc_track
        if (electron.row >= 0 && ahdc_track.row >= 0) {
            double delta_phi = fabs(electron.phi - ahdc_track.phi) - 180;
            H1_delta_phi_sym->Fill(delta_phi);
            H1_delta_phi_nosym->Fill(electron.phi - ahdc_track.phi);
            /////////////////////////////
            // elastics event
            // cut on delta_phi
            ////////////////////////////
            if (fabs(delta_phi) <= delta_phi_width) {
                H1_delta_phi_sym_sel->Fill(delta_phi);
                H1_delta_phi_nosym_sel->Fill(electron.phi - ahdc_track.phi);
                // correlations
                double pTe = 1000*sqrt(pow(electron.px,2) + pow(electron.py,2)); // MeV
                double pTt = sqrt(pow(ahdc_track.px,2) + pow(ahdc_track.py,2)); // MeV
                double sum_adc = trackBank.get("sum_adc", ahdc_track.row);
                H2_corr_pTe_Sadc->Fill(pTe,sum_adc);
                H2_corr_pTt_Sadc->Fill(pTt,sum_adc);
                // select the hits of this track
                int trackid = trackBank.getInt("trackid", ahdc_track.row);
                std::vector<int> hitRows;
                for (int hitRow = 0; hitRow < hitBank.getRows(); hitRow++) {
                    if (hitBank.getInt("trackid", hitRow) == trackid) {
                            hitRows.push_back(hitRow);
                    }
                }
                

                //////////////////
                // all
                //////////////////
                fill_histogram({H1_track_residual[0]}, {"residual"}, hitBank, hitRows);
                H2_corr_phi[0]->Fill(electron.phi, ahdc_track.phi);
                H2_corr_pT[0]->Fill(pTe, pTt); // Mev
                H2_corr_vz[0]->Fill(electron.vz, 0.1*ahdc_track.vz); // cm
                H1_delta_vz[0]->Fill(electron.vz-0.1*ahdc_track.vz); // cm
                H1_delta_phi[0]->Fill(delta_phi);
                // 1D
                H1_track_pT[0]->Fill(pTt);
                H1_track_phi[0]->Fill(ahdc_track.phi);
                H1_track_theta[0]->Fill(ahdc_track.theta);
                H1_electron_pT[0]->Fill(pTe);
                H1_electron_phi[0]->Fill(electron.phi);
                H1_electron_theta[0]->Fill(electron.theta);
                H1_diff_pT[0]->Fill(pTt - sqrt(pow(expected_track.px, 2) + pow(expected_track.py, 2)));
                H1_diff_theta[0]->Fill(ahdc_track.theta - expected_track.theta);
                H1_diff_phi[0]->Fill(ahdc_track.phi - expected_track.phi);
                //////////////////
                // deuteron cut
                //////////////////
                if (pTe > 200 && pTe < 300 && sum_adc > 5000 && sum_adc < 14000) {
                    fill_histogram({H1_track_residual[1]}, {"residual"}, hitBank, hitRows);
                    // 2D
                    H2_corr_phi[1]->Fill(electron.phi, ahdc_track.phi);
                    H2_corr_pT[1]->Fill(pTe, pTt); // Mev
                    H2_corr_vz[1]->Fill(electron.vz, 0.1*ahdc_track.vz); // cm
                    H1_delta_vz[1]->Fill(electron.vz-0.1*ahdc_track.vz); // cm
                    H1_delta_phi[1]->Fill(delta_phi);
                    // 1D
                    H1_track_pT[1]->Fill(pTt);
                    H1_track_phi[1]->Fill(ahdc_track.phi);
                    H1_track_theta[1]->Fill(ahdc_track.theta);
                    H1_electron_pT[1]->Fill(pTe);
                    H1_electron_phi[1]->Fill(electron.phi);
                    H1_electron_theta[1]->Fill(electron.theta);
                    H1_diff_pT[1]->Fill(pTt - sqrt(pow(expected_track.px, 2) + pow(expected_track.py, 2)));
                    H1_diff_theta[1]->Fill(ahdc_track.theta - expected_track.theta);
                    H1_diff_phi[1]->Fill(ahdc_track.phi - expected_track.phi);
                    
                }
                //////////////////
                // proton cut
                //////////////////
                if (pTe > 200 && pTe < 300 && sum_adc > 1000 && sum_adc < 4500) {
                    fill_histogram({H1_track_residual[2]}, {"residual"}, hitBank, hitRows);
                    // 2D
                    H2_corr_phi[2]->Fill(electron.phi, ahdc_track.phi);
                    H2_corr_pT[2]->Fill(pTe, pTt); // Mev
                    H2_corr_vz[2]->Fill(electron.vz, 0.1*ahdc_track.vz); // cm
                    H1_delta_vz[2]->Fill(electron.vz-0.1*ahdc_track.vz); // cm
                    H1_delta_phi[2]->Fill(delta_phi);
                    // 1D
                    H1_track_pT[2]->Fill(pTt);
                    H1_track_phi[2]->Fill(ahdc_track.phi);
                    H1_track_theta[2]->Fill(ahdc_track.theta);
                    H1_electron_pT[2]->Fill(pTe);
                    H1_electron_phi[2]->Fill(electron.phi);
                    H1_electron_theta[2]->Fill(electron.theta);
                    H1_diff_pT[2]->Fill(pTt - sqrt(pow(expected_track.px, 2) + pow(expected_track.py, 2)));
                    H1_diff_theta[2]->Fill(ahdc_track.theta - expected_track.theta);
                    H1_diff_phi[2]->Fill(ahdc_track.phi - expected_track.phi);   
                }
            }
        }

    }

    // save the output in a ROOT file
    TFile *f = new TFile(output.c_str(), "RECREATE");
    TCanvas* c_W2 = superimpose_histograms(H1_W2_all, H1_W2_sel, "canvas_W2");
    TCanvas* c_Q2 = superimpose_histograms(H1_Q2_all, H1_Q2_sel, "canvas_Q2");
    TCanvas* c_nu = superimpose_histograms(H1_nu_all, H1_nu_sel, "canvas_nu");
    TCanvas* c_delta_phi_sym = superimpose_histograms(H1_delta_phi_sym, H1_delta_phi_sym_sel, "canvas_delta_phi_sym");
    TCanvas* c_delta_phi_nosym = superimpose_histograms(H1_delta_phi_nosym, H1_delta_phi_nosym_sel, "canvas_delta_phi_nosym");
    TCanvas* c_elastics = display_elastics_cuts(H2_corr_pTe_Sadc, "canvas_elastics");
    c_W2->Write("W2");
    c_Q2->Write("Q2");
    c_nu->Write("nu");
    H1_vz->Write("electron_vz_nocuts");
    c_delta_phi_sym->Write("delta_phi_sym");
    c_delta_phi_nosym->Write("delta_phi_nosym");
    c_elastics->Write("corr_pT(electron)_sum_adc");
    H2_corr_pTt_Sadc->Write("corr_pT(track)_sum_adc");
    // all
    TDirectory *all_dir = f->mkdir("all_elastics");
    all_dir->cd();
    H2_corr_phi[0]->Write("corr_phi");
    H2_corr_pT[0]->Write("corr_pT");
    H2_corr_vz[0]->Write("corr_vz");
    TCanvas* c_delta_vz_0 = fit_histogram(H1_delta_vz[0], "canvas_delta_vz_0"); c_delta_vz_0->Write("delta_vz");
    TCanvas* c_delta_phi_0 = fit_histogram(H1_delta_phi[0], "canvas_delta_phi_0"); c_delta_phi_0->Write("delta_phi");
    H1_track_pT[0]->Write("h1_track_pT");
    H1_track_phi[0]->Write("h1_track_phi");
    H1_track_theta[0]->Write("h1_track_theta");
    H1_electron_pT[0]->Write("h1_electron_pT");
    H1_electron_phi[0]->Write("h1_electron_phi");
    H1_electron_theta[0]->Write("h1_electron_theta");
    H1_diff_pT[0]->Write("h1_diff_pT");
    H1_diff_theta[0]->Write("h1_diff_theta");
    H1_diff_phi[0]->Write("h1_diff_phi");
    //H1_track_residual[0]->Write("residual");
    TCanvas* c_residual_0 = fit_histogram(H1_track_residual[0], "canvas_residual_0"); c_residual_0->Write("residual");
    // deuteron
    TDirectory *deuteron_dir = f->mkdir("deuterons");
    deuteron_dir->cd();
    H2_corr_phi[1]->Write("corr_phi");
    H2_corr_pT[1]->Write("corr_pT");
    H2_corr_vz[1]->Write("corr_vz");
    TCanvas* c_delta_vz_1 = fit_histogram(H1_delta_vz[1], "canvas_delta_vz_1"); c_delta_vz_1->Write("delta_vz");
    TCanvas* c_delta_phi_1 = fit_histogram(H1_delta_phi[1], "canvas_delta_phi_1"); c_delta_phi_1->Write("delta_phi");
    H1_track_pT[1]->Write("h1_track_pT");
    H1_track_phi[1]->Write("h1_track_phi");
    H1_track_theta[1]->Write("h1_track_theta");
    H1_electron_pT[1]->Write("h1_electron_pT");
    H1_electron_phi[1]->Write("h1_electron_phi");
    H1_electron_theta[1]->Write("h1_electron_theta");
    H1_diff_pT[1]->Write("h1_diff_pT");
    H1_diff_theta[1]->Write("h1_diff_theta");
    H1_diff_phi[1]->Write("h1_diff_phi");
    //H1_track_residual[1]->Write("residual");
    TCanvas* c_residual_1 = fit_histogram(H1_track_residual[1], "canvas_residual_1"); c_residual_1->Write("residual");
    // proton
    TDirectory *proton_dir = f->mkdir("protons");
    proton_dir->cd();
    H2_corr_phi[2]->Write("corr_phi");
    H2_corr_pT[2]->Write("corr_pT");
    H2_corr_vz[2]->Write("corr_vz");
    TCanvas* c_delta_vz_2 = fit_histogram(H1_delta_vz[2], "canvas_delta_vz_2"); c_delta_vz_2->Write("delta_vz");
    TCanvas* c_delta_phi_2 = fit_histogram(H1_delta_phi[2], "canvas_delta_phi_2"); c_delta_phi_2->Write("delta_phi");
    H1_track_pT[2]->Write("h1_track_pT");
    H1_track_phi[2]->Write("h1_track_phi");
    H1_track_theta[2]->Write("h1_track_theta");
    H1_electron_pT[2]->Write("h1_electron_pT");
    H1_electron_phi[2]->Write("h1_electron_phi");
    H1_electron_theta[2]->Write("h1_electron_theta");
    H1_diff_pT[2]->Write("h1_diff_pT");
    H1_diff_theta[2]->Write("h1_diff_theta");
    H1_diff_phi[2]->Write("h1_diff_phi");
    //H1_track_residual[2]->Write("residual");
    TCanvas* c_residual_2 = fit_histogram(H1_track_residual[2], "canvas_residual_2"); c_residual_2->Write("residual");
    //H2_corr_pTe_Sadc->Write("corr_pT_Sadc");
    //H1_W2_all->Write("W2_all");
    //H1_delta_phi_nosym->Write("delta_phi_nosym_all");
    //H1_delta_phi_sym->Write("delta_phi_sym_all");

    // routine to create a directory in a TFile
    // TDirectory *adc_dir = f->mkdir("adc");
    // adc_dir->cd();
    // save a histogram
    //H1_mctime->Write("time");
    // ou H1_time->Write(H1_time->GetName());

    // small atof study
    count_atof_matching(filename, f);
    // close file
    f->Close();
    printf("File created : %s\n", output.c_str());
    
    // end of the program
    printf("* nevents : %ld\n", nevents);
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

TCanvas* superimpose_histograms(TH1D* h_all, TH1D* h_sel, const char * name) {
    TCanvas* c1 = new TCanvas(name);
    c1->cd();
    h_all->SetFillColor(kGray);
    h_sel->SetFillColor(kRed);
    h_all->Draw("hist");
    h_sel->Draw("hist same");
    return c1;
}

TCanvas* display_elastics_cuts(TH2D* h, const char * name) {
    TCanvas* c1 = new TCanvas(name);
    c1->cd();
    h->Draw("colz");
    { // deuterons
        double p1 = 200;
        double p2 = 300;
        double sum_adc1 = 5000;
        double sum_adc2 = 14000;
        TGraph * gr = new TGraph();
        gr->AddPoint(p1, sum_adc1);
        gr->AddPoint(p2, sum_adc1);
        gr->AddPoint(p2, sum_adc2);
        gr->AddPoint(p1, sum_adc2);
        gr->AddPoint(p1, sum_adc1);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(4);
        gr->Draw("same l");
    }
    { // protons
        double p1 = 200;
        double p2 = 300;
        double sum_adc1 = 1000;
        double sum_adc2 = 4500;
        TGraph * gr = new TGraph();
        gr->AddPoint(p1, sum_adc1);
        gr->AddPoint(p2, sum_adc1);
        gr->AddPoint(p2, sum_adc2);
        gr->AddPoint(p1, sum_adc2);
        gr->AddPoint(p1, sum_adc1);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(4);
        gr->Draw("same l");
    }
    return c1;
}

TCanvas* fit_histogram(TH1D* h, const char * name) {
    printf("> %s\n", name);
    TCanvas* c1 = new TCanvas(name);
    c1->cd();
    //h->SetFillColor(kGray);
    h->Draw("hist");
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    TF1 *f1 = new TF1("f1", "gaus(0) + gaus(3) + pol0(6)", xmin, xmax);
    //f1->SetLineColor(kGreen+2);
    f1->SetLineWidth(4);
    //f1->SetNpx(300);
    /*f1->SetParameter(0,h->GetBinContent(h->GetMaximumBin())); // amplitude
    f1->SetParameter(1,h->GetMean()); // mean 
    f1->SetParameter(2,h->GetStdDev()); // stdev
    f1->SetParameter(3,h->GetBinContent(h->GetMaximumBin())); // amplitude
    f1->SetParameter(4,h->GetMean()); // mean
    f1->SetParameter(5,h->GetStdDev()); // stdev
    f1->SetParameter(6,1); // baseline
    */
    
    f1->SetParameter(0,1); // amplitude
    f1->SetParameter(1,1); // mean 
    f1->SetParameter(2,1); // stdev
    f1->SetParameter(3,1); // amplitude
    f1->SetParameter(4,1); // mean
    f1->SetParameter(5,1); // stdev
    f1->SetParameter(6,1); // baseline
    
    h->Fit("f1", "", "", -xmin, xmax);
    gPad->Update();
    // get parameters
    double amp1 = f1->GetParameter(0);
    double mu1 = f1->GetParameter(1);
    double sigma1 = f1->GetParameter(2);
    double amp2 = f1->GetParameter(3);
    double mu2 = f1->GetParameter(4);
    double sigma2 = f1->GetParameter(5);
    //double baseline = f1->GetParameter(6);
    // compute stats
    //double mean = mu1 + mu2 + 0.5*(xmin+xmax);
    //double stdev = sqrt(pow(sigma1, 2) + pow(sigma2, 2));
    double mean = (amp1*mu1 + amp2*mu2)/(amp1+amp2);
    double mom1 = pow(sigma1,2) + pow(mu1,2); // moment of order 2: E(X^2)
    double mom2 = pow(sigma2,2) + pow(mu2,2);
    double mom = (amp1*mom1 + amp2*mom2)/(amp1+amp2);
    double stdev = sqrt(mom - pow(mean,2));
    double x1 = mean - stdev;
    double x2 = mean + stdev;
    // print info
    double half_amplitude = 0.5*f1->GetMaximum();
    //double x_for_max_value = f1->GetMaximumX();
    /*double half_amplitude = 0.5*(baseline + f1->GetMaximum());
    double x_for_max_value = f1->GetMaximumX();
    double x1 = f1->GetX(half_amplitude, -xmin, x_for_max_value);
    double x2 = f1->GetX(half_amplitude, x_for_max_value, xmax);*/

    TText* text = new TText();
    // SetTextAlign(10*h+v)
    // h horizontal:  1=left adjusted, 2=centered, 3=right adjusted
    // v vertical: 1=bottom adjusted, 2=centered, 3=top adjusted
    text->SetTextSize(0.03);
    text->SetTextAlign(12);
    text->SetTextColor(kGreen+2);
    text->DrawText(x2 + 0.05*(xmax-xmin), 1*half_amplitude, TString::Format("mean: %.3lf , error: %.3lf", 0.5*(x1+x2), 0.5*(x2-x1)).Data());
    // Draw width
    /*TGraph * gr = new TGraph();
    gr->AddPoint(x1, 0);
    gr->AddPoint(x1, half_amplitude);
    gr->AddPoint(0.5*(x1+x2), half_amplitude);
    gr->AddPoint(x2, half_amplitude);
    gr->AddPoint(x2, 0);
    gr->SetLineStyle(2);
    gr->SetLineColor(kRed);
    gr->Draw("same l");*/
    // Error graph
    TGraphErrors* gerr = new TGraphErrors();
    gerr->SetLineColor(kGreen+2);
    gerr->SetLineWidth(3);
    gerr->SetMarkerColor(kGreen+2);
    //gerr->SetMarkerSize(3);
    gerr->SetMarkerStyle(8);
    gerr->AddPointError(mean, 0.5*(f1->Eval(x1) + f1->Eval(x2)), stdev, 0);
    gerr->Draw("same lp");
    return c1;
}
void fill_histogram(std::vector<TH1D*> histos, std::vector<const char *> entries, hipo::bank & bank, std::vector<int> & rows) {
    // we should have the same number of histograms as for entries
    int nentries = histos.size();
    if (nentries != (int) entries.size()) {return;}
    for (int row : rows) {
        for (int e = 0; e < nentries; e++) {
            histos[e]->Fill(bank.get(entries[e], row));
        }
    }
}

void count_atof_matching(std::string filename, TFile * rootFile) {
    printf("////////////////////////////////////////\n");
    printf("/// Count ATOF matching\n");
    printf("////////////////////////////////////////\n");
    hipo::reader  reader(filename.c_str());
    hipo::dictionary factory;
    reader.readDictionary(factory);

    TH1D* H1_atof_bar_time = new TH1D("atof_bar_time", "time; time (ns); count", 50, 70, 112);
    TH1D* H1_atof_wedge_time = new TH1D("atof_wedge_time", "time; time (ns); count", 50, 70, 112);

    // bank definition
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  atofHitBank(factory.getSchema("ATOF::hits"));
    hipo::bank  aiMatchingBank(factory.getSchema("ALERT::ai:projections"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int ntracks =0;
    long unsigned int nmatches =0;

     while( reader.next()){
        nevents++;
        // display progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        // load bank content for this event
        reader.read(event);
        event.getStructure(trackBank);
        event.getStructure(atofHitBank);
        event.getStructure(aiMatchingBank);
        ///////////////////////////////////
        /// count ahdc aof matches
        ///////////////////////////////////
        if (trackBank.getRows() > 0 && atofHitBank.getRows() > 0) {
            //ntracks += trackBank.getRows();
            ntracks++;
            bool isFirst = false;
            for (int aiRow = 0; aiRow < aiMatchingBank.getRows(); aiRow++) {
                int trackid = aiMatchingBank.get("trackid", aiRow);
                int atofid = aiMatchingBank.get("matched_atof_hit_id", aiRow);
                if (trackid > 0 && atofid > 0 && !isFirst) {
                    nmatches++;
                    isFirst = true;
                    //aiMatchingBank.show();
                }
            }
        }
        //////////////////////////////////////
        /// atof bar time
        //////////////////////////////////////
        for (int hitRow = 0; hitRow < atofHitBank.getRows(); hitRow++) {
            int component = atofHitBank.get("component", hitRow);
            // bar selection
            double time = atofHitBank.get("time", hitRow);
            if (component == 10) {
                H1_atof_bar_time->Fill(time);
            } // else are wedges
            else {
                H1_atof_wedge_time->Fill(time);
            }
        }

    }

    TDirectory *atof_dir = rootFile->mkdir("atof");
    atof_dir->cd();
    H1_atof_bar_time->Write("bar_time");
    H1_atof_wedge_time->Write("wedge_time");
    printf("      Number of events : %ld\n", nevents);
    printf("    Number of kftracks : %ld   ===>  %lf %%\n", ntracks, 100.0*ntracks/nevents);
    printf("Number of atof matches : %ld   ===>  %lf %%\n", nmatches, 100.0*nmatches/ntracks);
}

#include <filesystem>
bool read_options(int argc, const char * argv, std::string filename, std::string outname) {
    if (argc < 2) {
        printf("Please provide 1 file name and 1 output name\n");
        return false;
    }
    /*if (std::to_string(argv[1]).compare("-f") == 0) {
        // require that the next argument should be a filename
        filename = argv[2];
        if (!std::filesystem::is_regular_file(filename.c_str())) {
            printf("%s does not exist.\n", filename.c_str());
            return false;
        }
        else {
            if ()
        }

    }*/
   return true;
}