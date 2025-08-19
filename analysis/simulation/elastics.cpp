/***********************************************
 * Analysis of tracks	
 *
 * @author Felix Touchte Codjo
 * @date May 21, 2025
 * ********************************************/

/** Notes on the cuts
 * 0) requires at least one track, startTime > 0, one electron 
 * Exclude 
 * 1) electron -> vz < -25 and vz > 10
 * 2) track -> theta < 10 and theta > 170 
 *      - no number of hits == 2
 *      - always nhits from 4
 *      - always path at 0
 *      - vz for track always above 15 cm
 * 3) track -> |vz| > 15 cm
 *      - always nhits from 4
 *      - almost no longer path at 0
 * 4) track -> nhits >= 6
 * 5) W2 < 14 and W2 > 13.85
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>
#include <filesystem>

#include "../hipo4/reader.h"
#include "../hipo4/writer.h"

#include "futils.h"
#include "elastics.h"
#include "AhdcCCDB.h"

#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "THStack.h"


// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

const double Ee = 2.23951; // incident energy if the electron, GeV
const double me = 0.511e-3; // energy mass of electron, GeV
const double Mp = 938.272e-3; // energy mass of proton, GeV
const double M_He = 3.73; // energy mass of Helium-4, GeV
const double M_D = 1.875; // energy mass of Deuterium, GeV
const double Mt = M_D; // target rest mass : choose Mp or M_He, GeV

// Physics - Histogram limits
const double lim_W2_inf = 2;
const double lim_W2_sup = 7;
const double lim_Q2_inf = 0;
const double lim_Q2_sup = 0.4;
const double lim_xB_inf = 0;
const double lim_xB_sup = 2;
const double lim_nu_inf = 0;
const double lim_nu_sup = 1.9;
// Probe - Electron or Photon
const double lim_probe_p_inf = 0;
const double lim_probe_p_sup = 2.5;
const double lim_probe_pT_inf = 0;
const double lim_probe_pT_sup = 0.9;
const double lim_probe_theta_inf = 2;
const double lim_probe_theta_sup = 25;
const double lim_probe_phi_inf = 0;
const double lim_probe_phi_sup = 361;
const double lim_probe_vz_inf = -45;
const double lim_probe_vz_sup = 35;
// AHDC KF track 
const double lim_track_p_inf = 0;
const double lim_track_p_sup = 1.5;
const double lim_track_pT_inf = 0.2;
const double lim_track_pT_sup = 0.8;
const double lim_track_theta_inf = 0;
const double lim_track_theta_sup = 181;
const double lim_track_phi_inf = 0;
const double lim_track_phi_sup = 361;
const double lim_track_vz_inf = -30;
const double lim_track_vz_sup = 30;


int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    double W2_min = 3.46;
    double W2_max = 3.67;
    //double delta_W2 = W2_max - W2_min;
    //////////////////////////////
    // Open HIPO file
    //////////////////////////////
    // He4
    //const char * filename = "/home/touchte-codjo/Desktop/hipofiles/track/He4/half_det_lost/all_half_det_lost.hipo";
    //const char * filename = "/home/touchte-codjo/Desktop/hipofiles/track/He4/22767/all_rec_clas_022767.hipo";
    // D2
    //const char * filename = "/home/touchte-codjo/Desktop/hipofiles/track/D2/22712/rec_clas_022712.evio.00000.hipo"; 
    const char * filename = "/home/touchte-codjo/Desktop/hipofiles/track/D2/22712/all_rec_clas_022712.hipo";

    AhdcCCDB ahdcConstants; // load ccdb
    // Input File
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  track0Bank(factory.getSchema("AHDC::track")); // to match old cooked files where dEdx was filled in AHDC::track instead of AHDC::kftrack
    hipo::bank  particleBank(factory.getSchema("REC::Particle"));
    hipo::bank  recEventBank(factory.getSchema("REC::Event"));
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::event event;
    // Ouput file
    hipo::writer writer;
    //factory.getSchema("AHDC::adc").show(); 
    //factory.getSchema("AHDC::wf").show(); 
    //factory.getSchema("AHDC::kftrack").show();
    hipo::schema schemaAdc("AHDC::adc", 22400, 11);
    hipo::schema schemaWf("AHDC::wf", 22400, 10);
    hipo::schema schemaTrack("AHDC::kftrack",23000, 26);
    schemaAdc.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/F, windex/S, integral/I, leadingEdgeTime/F, timeOverThreshold/F, constantFractionTime/F, wfType/S");
    schemaWf.parse("sector/B, layer/B, component/S, order/B, timestamp/L, s1/S, s2/S, s3/S, s4/S, s5/S, s6/S, s7/S, s8/S, s9/S, s10/S, s11/S, s12/S, s13/S, s14/S, s15/S, s16/S, s17/S, s18/S, s19/S, s20/S, s21/S, s22/S, s23/S, s24/S, s25/S, s26/S, s27/S, s28/S, s29/S, s30/S, time/I");
    schemaTrack.parse("trackid/I, x/F, y/F, z/F, px/F, py/F, pz/F, n_hits/I, sum_adc/I, path/F, dEdx/F, p_drift/F, chi2/F, sum_residuals/F");
    writer.getDictionary().addSchema(schemaAdc);
    writer.getDictionary().addSchema(schemaWf);
    writer.getDictionary().addSchema(schemaTrack);
    const char * outputFile = "elastics_events.hipo";
    if (std::filesystem::exists(outputFile)) {
        if (std::filesystem::remove(outputFile)) {
            printf("Remove file before processing : %s\n", outputFile);
        }
    }
    writer.open(outputFile);
    hipo::event outEvent;

    // Histograms
    long unsigned int nevents = 0;
    long unsigned int nelectrons = 0;
    long unsigned int nphotons = 0;
    long unsigned int ntracks = 0;
    TH1I* H1_mon_nprobes    = new TH1I("nprobes_with_cuts", "number of probes with the cuts", 15, 0, 15); 
    TH1I* H1_mon_ntracks    = new TH1I("ntracks_with_cuts", "number of tracks with the cuts", 15, 0, 15); 
    //////////////////////////////
    // Probe : electron or photon
    //////////////////////////////
    std::vector<TH1D*> H1_probe_vz;
    H1_probe_vz.push_back(new TH1D("vz_probe_0", "vz electron/photon (cm)", 100, lim_probe_vz_inf, lim_probe_vz_sup));
    H1_probe_vz.push_back(new TH1D("vz_probe_1", "vz electron/photon (cm)", 100, lim_probe_vz_inf, lim_probe_vz_sup));
    H1_probe_vz.push_back(new TH1D("vz_probe_2", "vz electron/photon (cm)", 100, lim_probe_vz_inf, lim_probe_vz_sup));
    H1_probe_vz.push_back(new TH1D("vz_probe_3", "vz electron/photon (cm)", 100, lim_probe_vz_inf, lim_probe_vz_sup));
    std::vector<TH1D*> H1_probe_p;
    H1_probe_p.push_back(new TH1D("p_probe_0", "p electron/photon (GeV)", 100, lim_probe_p_inf, lim_probe_p_sup));
    H1_probe_p.push_back(new TH1D("p_probe_1", "p electron/photon (GeV)", 100, lim_probe_p_inf, lim_probe_p_sup));
    H1_probe_p.push_back(new TH1D("p_probe_2", "p electron/photon (GeV)", 100, lim_probe_p_inf, lim_probe_p_sup));
    H1_probe_p.push_back(new TH1D("p_probe_3", "p electron/photon (GeV)", 100, 2.16, 2.25));
    std::vector<TH1D*> H1_probe_pT;
    H1_probe_pT.push_back(new TH1D("pT_probe_0", "pT electron/photon (GeV)", 100, lim_probe_pT_inf, lim_probe_pT_sup));
    H1_probe_pT.push_back(new TH1D("pT_probe_1", "pT electron/photon (GeV)", 100, lim_probe_pT_inf, lim_probe_pT_sup));
    H1_probe_pT.push_back(new TH1D("pT_probe_2", "pT electron/photon (GeV)", 100, lim_probe_pT_inf, lim_probe_pT_sup));
    H1_probe_pT.push_back(new TH1D("pT_probe_3", "pT electron/photon (GeV)", 100, 0.1, 0.5));
    std::vector<TH1D*> H1_probe_theta;
    H1_probe_theta.push_back(new TH1D("theta_probe_0", "theta electron/photon (deg)", 100, lim_probe_theta_inf, lim_probe_theta_sup));
    H1_probe_theta.push_back(new TH1D("theta_probe_1", "theta electron/photon (deg)", 100, lim_probe_theta_inf, lim_probe_theta_sup));
    H1_probe_theta.push_back(new TH1D("theta_probe_2", "theta electron/photon (deg)", 100, lim_probe_theta_inf, lim_probe_theta_sup));
    H1_probe_theta.push_back(new TH1D("theta_probe_3", "theta electron/photon (deg)", 100, 2.56, 12));
    std::vector<TH1D*> H1_probe_phi;
    H1_probe_phi.push_back(new TH1D("phi_probe_0", "phi electron/photon (deg)", 100, lim_probe_phi_inf, lim_probe_phi_sup));
    H1_probe_phi.push_back(new TH1D("phi_probe_1", "phi electron/photon (deg)", 100, lim_probe_phi_inf, lim_probe_phi_sup));
    H1_probe_phi.push_back(new TH1D("phi_probe_2", "phi electron/photon (deg)", 100, lim_probe_phi_inf, lim_probe_phi_sup));
    H1_probe_phi.push_back(new TH1D("phi_probe_3", "phi electron/photon (deg)", 100, lim_probe_phi_inf, lim_probe_phi_sup));
    std::vector<TH1I*> H1_nelectrons;
    H1_nelectrons.push_back(new TH1I("nelectrons_0", "number of electrons / evt", 5, 0, 5));
    H1_nelectrons.push_back(new TH1I("nelectrons_1", "number of electrons / evt", 5, 0, 5));
    H1_nelectrons.push_back(new TH1I("nelectrons_2", "number of electrons / evt", 5, 0, 5));
    H1_nelectrons.push_back(new TH1I("nelectrons_3", "number of electrons / evt", 5, 0, 5));
    std::vector<TH1I*> H1_nphotons;
    H1_nphotons.push_back(new TH1I("nphotons_0", "number of photons / evt", 5, 0, 5));
    H1_nphotons.push_back(new TH1I("nphotons_1", "number of photons / evt", 5, 0, 5));
    H1_nphotons.push_back(new TH1I("nphotons_2", "number of photons / evt", 5, 0, 5));
    H1_nphotons.push_back(new TH1I("nphotons_3", "number of photons / evt", 5, 0, 5));
    //////////////////////////////
    // physics
    //////////////////////////////
    std::vector<TH1D*> H1_Q2;
    H1_Q2.push_back(new TH1D("Q2_0", "Q^{2} (GeV^{2})", 100, lim_Q2_inf, lim_Q2_sup));
    H1_Q2.push_back(new TH1D("Q2_1", "Q^{2} (GeV^{2})", 100, lim_Q2_inf, lim_Q2_sup));
    H1_Q2.push_back(new TH1D("Q2_2", "Q^{2} (GeV^{2})", 100, lim_Q2_inf, lim_Q2_sup));
    H1_Q2.push_back(new TH1D("Q2_3", "Q^{2} (GeV^{2})", 100, lim_Q2_inf, lim_Q2_sup));
    std::vector<TH1D*> H1_W2;
    H1_W2.push_back(new TH1D("W2_0", "W^{2} (GeV^{2})", 100, lim_W2_inf, lim_W2_sup));
    H1_W2.push_back(new TH1D("W2_1", "W^{2} (GeV^{2})", 100, lim_W2_inf, lim_W2_sup));
    H1_W2.push_back(new TH1D("W2_2", "W^{2} (GeV^{2})", 100, lim_W2_inf, lim_W2_sup));
    H1_W2.push_back(new TH1D("W2_3", "W^{2} (GeV^{2})", 100, lim_W2_inf, lim_W2_sup));
    std::vector<TH1D*> H1_xB;
    H1_xB.push_back(new TH1D("xB_0", "x_{B}", 100, lim_xB_inf, lim_xB_sup));
    H1_xB.push_back(new TH1D("xB_1", "x_{B}", 100, lim_xB_inf, lim_xB_sup));
    H1_xB.push_back(new TH1D("xB_2", "x_{B}", 100, lim_xB_inf, lim_xB_sup));
    H1_xB.push_back(new TH1D("xB_3", "x_{B}", 100, 0.15, 4));
    std::vector<TH1D*> H1_nu;
    H1_nu.push_back(new TH1D("nu_0", "#nu = #Delta E = E - E' (GeV)", 100, lim_nu_inf, lim_nu_sup)); 
    H1_nu.push_back(new TH1D("nu_1", "#nu = #Delta E = E - E' (GeV)", 100, lim_nu_inf, lim_nu_sup)); 
    H1_nu.push_back(new TH1D("nu_2", "#nu = #Delta E = E - E' (GeV)", 100, lim_nu_inf, lim_nu_sup)); 
    H1_nu.push_back(new TH1D("nu_3", "#nu = #Delta E = E - E' (GeV)", 100, lim_nu_inf, 0.1)); 
    //////////////////////////////
    // ahdc (kf) track
    //////////////////////////////
    std::vector<TH1I*> H1_ntracks;
    H1_ntracks.push_back(new TH1I("ntracks_0", "number of tracks / evt", 5, 0, 5));
    H1_ntracks.push_back(new TH1I("ntracks_1", "number of tracks / evt", 5, 0, 5));
    H1_ntracks.push_back(new TH1I("ntracks_2", "number of tracks / evt", 5, 0, 5));
    H1_ntracks.push_back(new TH1I("ntracks_3", "number of tracks / evt", 5, 0, 5));
    std::vector<TH1D*> H1_track_vz;
    H1_track_vz.push_back(new TH1D("vz_track_0", "vz (cm)", 100, -30, 30));
    H1_track_vz.push_back(new TH1D("vz_track_1", "vz (cm)", 100, -30, 30));
    H1_track_vz.push_back(new TH1D("vz_track_2", "vz (cm)", 100, -30, 30));
    H1_track_vz.push_back(new TH1D("vz_track_3", "vz (cm)", 100, -30, 30));
    std::vector<TH1D*> H1_track_p;
    H1_track_p.push_back(new TH1D("p_track_0", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p.push_back(new TH1D("p_track_1", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p.push_back(new TH1D("p_track_2", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p.push_back(new TH1D("p_track_3", "p (GeV)", 100, 0.2, 1.6));
    std::vector<TH1D*> H1_track_p_drift;
    H1_track_p_drift.push_back(new TH1D("p_drift_track_0", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p_drift.push_back(new TH1D("p_drift_track_1", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p_drift.push_back(new TH1D("p_drift_track_2", "p (GeV)", 100, 0.2, 1.6));
    H1_track_p_drift.push_back(new TH1D("p_drift_track_3", "p (GeV)", 100, 0.2, 1.6));
    std::vector<TH1D*> H1_track_pT;
    H1_track_pT.push_back(new TH1D("pT_track_0", "pT (GeV)", 100, 0, 1));
    H1_track_pT.push_back(new TH1D("pT_track_1", "pT (GeV)", 100, 0, 1));
    H1_track_pT.push_back(new TH1D("pT_track_2", "pT (GeV)", 100, 0, 1));
    H1_track_pT.push_back(new TH1D("pT_track_3", "pT (GeV)", 100, 0, 1));
    std::vector<TH1D*> H1_track_theta;
    H1_track_theta.push_back(new TH1D("theta_track_0", "theta (deg)", 100, 0, 181));
    H1_track_theta.push_back(new TH1D("theta_track_1", "theta (deg)", 100, 0, 181));
    H1_track_theta.push_back(new TH1D("theta_track_2", "theta (deg)", 100, 0, 181));
    H1_track_theta.push_back(new TH1D("theta_track_3", "theta (deg)", 100, 0, 181));
    std::vector<TH1D*> H1_track_phi;
    H1_track_phi.push_back(new TH1D("phi_track_0", "phi (deg)", 100, 0, 361));
    H1_track_phi.push_back(new TH1D("phi_track_1", "phi (deg)", 100, 0, 361));
    H1_track_phi.push_back(new TH1D("phi_track_2", "phi (deg)", 100, 0, 361));
    H1_track_phi.push_back(new TH1D("phi_track_3", "phi (deg)", 100, 0, 361));
    std::vector<TH1D*> H1_track_nhits;
    H1_track_nhits.push_back(new TH1D("nhits_track_0", "nhits per track (deg)", 100, 0, 15));
    H1_track_nhits.push_back(new TH1D("nhits_track_1", "nhits per track (deg)", 100, 0, 15));
    H1_track_nhits.push_back(new TH1D("nhits_track_2", "nhits per track (deg)", 100, 0, 15));
    H1_track_nhits.push_back(new TH1D("nhits_track_3", "nhits per track (deg)", 100, 0, 15));
    std::vector<TH1D*> H1_track_adc;
    H1_track_adc.push_back(new TH1D("adc_track_0", "Sum adc (adc)", 100, 0, 25000));
    H1_track_adc.push_back(new TH1D("adc_track_1", "Sum adc (adc)", 100, 0, 25000));
    H1_track_adc.push_back(new TH1D("adc_track_2", "Sum adc (adc)", 100, 0, 25000));
    H1_track_adc.push_back(new TH1D("adc_track_3", "Sum adc (adc)", 100, 0, 25000));
    std::vector<TH1D*> H1_track_path;
    H1_track_path.push_back(new TH1D("path_track_0", "path (mm)", 100, 40, 300));
    H1_track_path.push_back(new TH1D("path_track_1", "path (mm)", 100, 40, 300));
    H1_track_path.push_back(new TH1D("path_track_2", "path (mm)", 100, 40, 300));
    H1_track_path.push_back(new TH1D("path_track_3", "path (mm)", 100, 40, 300));
    std::vector<TH1D*> H1_track_dEdx;
    H1_track_dEdx.push_back(new TH1D("dEdx_track_0", "dEdx (adc/mm)", 100, 0, 400));
    H1_track_dEdx.push_back(new TH1D("dEdx_track_1", "dEdx (adc/mm)", 100, 0, 400));
    H1_track_dEdx.push_back(new TH1D("dEdx_track_2", "dEdx (adc/mm)", 100, 0, 400));
    H1_track_dEdx.push_back(new TH1D("dEdx_track_3", "dEdx (adc/mm)", 100, 0, 400));
    std::vector<TH1D*> H1_track_residuals;
    H1_track_residuals.push_back(new TH1D("residuals_track_0", "Sum residuals(mm)", 100, -25, 3));
    H1_track_residuals.push_back(new TH1D("residuals_track_1", "Sum residuals(mm)", 100, -25, 3));
    H1_track_residuals.push_back(new TH1D("residuals_track_2", "Sum residuals(mm)", 100, -25, 3));
    H1_track_residuals.push_back(new TH1D("residuals_track_3", "Sum residuals(mm)", 100, -25, 3));
    std::vector<TH1D*> H1_track_residuals_per_nhits;
    H1_track_residuals_per_nhits.push_back(new TH1D("residuals_per_nhits_track_0", "residual per hit (mm)", 100, -4.17, 0.5));
    H1_track_residuals_per_nhits.push_back(new TH1D("residuals_per_nhits_track_1", "residual per hit (mm)", 100, -4.17, 0.5));
    H1_track_residuals_per_nhits.push_back(new TH1D("residuals_per_nhits_track_2", "residual per hit (mm)", 100, -4.17, 0.5));
    H1_track_residuals_per_nhits.push_back(new TH1D("residuals_per_nhits_track_3", "residual per hit (mm)", 100, -4.17, 0.5));
    std::vector<TH1D*> H1_track_chi2;
    H1_track_chi2.push_back(new TH1D("chi2_track_0", "chi2 (mm^{2})", 100, 0, 100));
    H1_track_chi2.push_back(new TH1D("chi2_track_1", "chi2 (mm^{2})", 100, 0, 100));
    H1_track_chi2.push_back(new TH1D("chi2_track_2", "chi2 (mm^{2})", 100, 0, 100));
    H1_track_chi2.push_back(new TH1D("chi2_track_3", "chi2 (mm^{2})", 100, 0, 100));
    std::vector<TH1D*> H1_track_chi2_per_nhits;
    H1_track_chi2_per_nhits.push_back(new TH1D("chi2_per_nhits_track_0", "chi2 per hit (mm^{2})", 100, 0, 16.67));
    H1_track_chi2_per_nhits.push_back(new TH1D("chi2_per_nhits_track_1", "chi2 per hit (mm^{2})", 100, 0, 16.67));
    H1_track_chi2_per_nhits.push_back(new TH1D("chi2_per_nhits_track_2", "chi2 per hit (mm^{2})", 100, 0, 16.67));
    H1_track_chi2_per_nhits.push_back(new TH1D("chi2_per_nhits_track_3", "chi2 per hit (mm^{2})", 100, 0, 16.67));
    //////////////////////////////////////////////////////
    // correlations probes (electron or gamma) vs tracks)
    //////////////////////////////////////////////////////
    std::vector<TH2D*> H2_p_dEdx;
    H2_p_dEdx.push_back(new TH2D("pTe_dEdx_0", "pT electron vs dEdx", 100, 0.2, 0.45, 100, 0, 200));
    H2_p_dEdx.push_back(new TH2D("pTe_dEdx_1", "pT electron vs dEdx", 100, 0.2, 0.45, 100, 0, 200));
    H2_p_dEdx.push_back(new TH2D("pTe_dEdx_2", "pT electron vs dEdx", 100, 0.2, 0.45, 100, 0, 200));
    H2_p_dEdx.push_back(new TH2D("pTe_dEdx_3", "pT electron vs dEdx", 100, 0.2, 0.45, 100, 0, 200));
    std::vector<TH2D*> H2_p_adc;
    H2_p_adc.push_back(new TH2D("pTe_adc_0", "pT e^{-}/#gamma vs #Sigma adc", 100, 0.19, 0.45, 100, 0, 15000));
    H2_p_adc.push_back(new TH2D("pTe_adc_1", "pT e^{-}/#gamma vs #Sigma adc", 100, 0.19, 0.45, 100, 0, 15000));
    H2_p_adc.push_back(new TH2D("pTe_adc_2", "pT e^{-}/#gamma vs #Sigma adc", 100, 0.19, 0.45, 100, 0, 15000));
    H2_p_adc.push_back(new TH2D("pTe_adc_3", "pT e^{-}/#gamma vs #Sigma adc", 100, 0.19, 0.45, 100, 0, 15000));
    std::vector<TH2D*> H2_vze_vz;
    H2_vze_vz.push_back(new TH2D("vze_vz_0", "vz_{e} vs vz", 100, -25, 10, 100, -16, 16)); 
    H2_vze_vz.push_back(new TH2D("vze_vz_1", "vz_{e} vs vz", 100, -25, 10, 100, -16, 16)); 
    H2_vze_vz.push_back(new TH2D("vze_vz_2", "vz_{e} vs vz", 100, -25, 10, 100, -16, 16)); 
    H2_vze_vz.push_back(new TH2D("vze_vz_3", "vz_{e} vs vz", 100, -25, 10, 100, -16, 16)); 
    std::vector<TH2D*> H2_pTe_pT;
    H2_pTe_pT.push_back(new TH2D("pTe_pT_0", "pT_{e} vs pT", 100, 0, 0.9, 100, 0, 1)); 
    H2_pTe_pT.push_back(new TH2D("pTe_pT_1", "pT_{e} vs pT", 100, 0, 0.9, 100, 0, 1)); 
    H2_pTe_pT.push_back(new TH2D("pTe_pT_2", "pT_{e} vs pT", 100, 0, 0.9, 100, 0, 1)); 
    H2_pTe_pT.push_back(new TH2D("pTe_pT_3", "pT_{e} vs pT", 100, 0, 0.9, 100, 0, 1)); 
    std::vector<TH1D*> H1_delta_vz;
    H1_delta_vz.push_back(new TH1D("delta_vz_0", "#Delta vz = vz_{e} - vz (cm)", 100, -40, 40));
    H1_delta_vz.push_back(new TH1D("delta_vz_1", "#Delta vz = vz_{e} - vz (cm)", 100, -40, 40));
    H1_delta_vz.push_back(new TH1D("delta_vz_2", "#Delta vz = vz_{e} - vz (cm)", 100, -40, 40));
    H1_delta_vz.push_back(new TH1D("delta_vz_3", "#Delta vz = vz_{e} - vz (cm)", 100, -40, 40));
    std::vector<TH1D*> H1_delta_phi; 
    H1_delta_phi.push_back(new TH1D("delta_phi_0", "#Delta #phi = #phi_{e} - #phi (deg)", 100, -360, 360));
    H1_delta_phi.push_back(new TH1D("delta_phi_1", "#Delta #phi = #phi_{e} - #phi (deg)", 100, -360, 360));
    H1_delta_phi.push_back(new TH1D("delta_phi_2", "#Delta #phi = #phi_{e} - #phi (deg)", 100, -360, 360));
    H1_delta_phi.push_back(new TH1D("delta_phi_3", "#Delta #phi = #phi_{e} - #phi (deg)", 100, -360, 360));
    //////////////////////////////////////////////////////
    /// for simulation
    //////////////////////////////////////////////////////
    TH1D* H1_t0 = new TH1D("t0", "t0; time (ns); count", 100, 150, 400); 
    TH1D* H1_dt0 = new TH1D("dt0", "#Delta t_0; leadingEdgeTime - t0 (ns); count", 100, 0, 400); 
    TH1D* H1_dt00 = new TH1D("dt0_before_cuts", "#Delta t_0 before cuts; leadingEdgeTime - t0 (ns); count", 100, 0, 400); 
    TH1I* H1_nelastics = new TH1I("nelastics", "nelastics per event; nelastics; count;", 10, 0, 10); 
    TH1D* H1_leadingEdgeTime = new TH1D("leadingEdgeTime", "leadingEdgeTime; time (ns); count", 100, 0, 700); 
    TH1D* H1_timeMax = new TH1D("timeMax", "timeMax; timeMax (ns); count", 100, 200, 900); 
    TH1D* H1_deltaTime = new TH1D("deltaTime", "timeMax - leadingEdgeTime (ns); timeMax - leadingEdgeTime (ns); count", 100, 0, 400); 
    TH1D* H1_timeOverThreshold = new TH1D("timeOverThreshold", "timeOverThreshol (ns); timeOverThreshol (ns); count", 100, 150, 752); 
    TH1D* H1_tot0 = new TH1D("timeOverThreshold0", "ToT just after dt0 cut; timeOverThreshol (ns); count", 100, 150, 752); 
    TH1I* H1_wfType = new TH1I("wfType", "wfType; count;", 6, 0, 6); 
    TH1I* H1_amplitude = new TH1I("amplitude", "amplitude (adc); count;", 100, 0, 3700); 
    TH2D* H2_times = new TH2D("timeMax, leadingEdgeTime", "timeMax vs leadingEdgeTime; timeMax (ns); leadingEdgeTime (ns);", 10, 200, 900, 100, 0, 700); 
    TH2D* H2_amp_tot = new TH2D("amp, tot", "amplitude vs timeOverThreshold;amp (adc); tot (ns)", 100, 0, 3700, 100, 150, 752); 
    TH2D* H2_deltaTime_adc = new TH2D("deltaTime_adc", "deltaTime vs amplitude;deltaTime (ns); amplitude (ns)", 100, 0, 400, 100, 0, 3700); 
    TH1D* H1_elastic_track_p = new TH1D("elastic_track_p", "p (GeV)", 100, 0.2, 1.6); 
    TH1D* H1_elastic_track_pT = new TH1D("elastic_track_pT", "pT (GeV)", 100, 0, 1); 
    TH1D* H1_elastic_track_theta = new TH1D("elastic_track_theta", "theta (deg)", 100, 0, 181); 
    TH1D* H1_elastic_track_phi = new TH1D("elastic_track_phi", "phi (deg)", 100, 0, 361); 
    TH1D* H1_elastic_probe_p = new TH1D("elastic_probe_p", "p (GeV)", 100, 2.16, 2.25); 
    TH1D* H1_elastic_probe_pT = new TH1D("elastic_probe_pT", "pT (GeV)", 100, 0.1, 0.5); 
    TH1D* H1_elastic_probe_theta = new TH1D("elastic_probe_theta", "theta (deg)", 100, 0, 12); 
    TH1D* H1_elastic_probe_phi = new TH1D("elastic_probe_phi", "phi (deg)", 100, 0, 361);
    TH2D* H2_elastic_pT_adc = new TH2D("elastic_pT_adc", "pT e^{-}/#gamma vs #Sigma adc; pT_{probe}; #Sigma adc", 100, 0.19, 0.45, 100, 0, 15000);
    TH1D* H1_elastic_theory_p = new TH1D("theory_track_p", "p (GeV)", 100, 0.2, 1.6); // what should be the track knowing it is an elastic and from theta_electron 
    TH1D* H1_elastic_theory_pT = new TH1D("theory_track_pT", "pT (GeV)", 100, 0, 1); 
    TH1D* H1_elastic_theory_theta = new TH1D("theory_track_theta", "theta (deg)", 100, 0, 181); 
    TH1D* H1_elastic_theory_phi = new TH1D("theory_track_phi", "phi (deg)", 100, 0, 361); 
    TH1D* H1_selection_theory_p = new TH1D("selection_theory_track_p", "p (GeV)", 100, 0.229, 0.271); 
    TH1D* H1_selection_theory_pT = new TH1D("selection_theory_track_pT", "pT (GeV)", 100, 0.22, 0.27); 
    TH1D* H1_selection_theory_theta = new TH1D("selection_theory_track_theta", "theta (deg)", 100, 86.55, 87.05); 
    TH1D* H1_selection_theory_phi = new TH1D("selection_theory_track_phi", "phi (deg)", 100, 0, 361); 
    TH1D* H1_selection_theory_adc = new TH1D("selection_theory_track_adc", "#Sum adc; #Sum adc; count", 100, 950, 3550); 
    // photons in FT
    TH1D* H1_photon_p = new TH1D("ft_photon_p", "FT/photon p (GeV); p (GeV); count", 100, 0, 2.5); 
    TH1D* H1_photon_pT = new TH1D("ft_photon_pT", "FT/photon pT (GeV); pT (GeV); count", 100, 0, 0.21); 
    TH1D* H1_photon_theta = new TH1D("ft_photon_theta", "FT/photon theta (deg); theta (deg); count", 100, 2, 5); 
    TH1D* H1_photon_phi = new TH1D("ft_photon_phi", "FT/photon phi (deg); phi (deg); count", 100, 0, 361); 

    // Loop over events
    while( reader.next()){
        nevents++;
        // Progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        reader.read(event);
        event.getStructure(trackBank);
        event.getStructure(track0Bank);
        event.getStructure(particleBank);
        event.getStructure(recEventBank);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(hitBank);
        
        //////////////////////
        //  REC::Particle
        //////////////////////
        std::vector<State> Electrons; 
        std::vector<State> Photons; 
        std::vector<Physics> ElectronKinematics;
        std::vector<Physics> PhotonKinematics;
        for (int i = 0; i < particleBank.getRows(); i++) {
            int pid = particleBank.getInt("pid", i);
            if (pid == 11) { // electron (11)
                State S(particleBank.getFloat("px",i), particleBank.getFloat("py",i), particleBank.getFloat("pz",i), 
                        particleBank.getFloat("vx",i), particleBank.getFloat("vy",i), particleBank.getFloat("vz",i));
                S.pid = 11;
                Physics K(S, Mt, Ee);
                Electrons.push_back(S);
                ElectronKinematics.push_back(K);
                // Histograms (0) no cuts
                H1_probe_vz[0]->Fill(S.vz);
                H1_probe_p[0]->Fill(S.p);
                H1_probe_pT[0]->Fill(S.pT);
                H1_probe_theta[0]->Fill(S.theta);
                H1_probe_phi[0]->Fill(S.phi);
                H1_Q2[0]->Fill(K.Q2);
                H1_W2[0]->Fill(K.W2);
                H1_xB[0]->Fill(K.xB);
                H1_nu[0]->Fill(K.nu);
            }
            else if (pid == 22) { // photon (22)
                State S(particleBank.getFloat("px",i), particleBank.getFloat("py",i), particleBank.getFloat("pz",i), 
                        particleBank.getFloat("vx",i), particleBank.getFloat("vy",i), particleBank.getFloat("vz",i));
                S.pid = 22;
                if (S.theta < 5) { // photon distribution in the FT
                    H1_photon_p->Fill(S.p);
                    H1_photon_pT->Fill(S.pT);
                    H1_photon_theta->Fill(S.theta);
                    H1_photon_phi->Fill(S.phi);
                }
                Physics K(S, Ee, Mt);
                Photons.push_back(S);
                PhotonKinematics.push_back(K);
                // Histograms (0) no cuts
                H1_probe_vz[0]->Fill(S.vz);
                H1_probe_p[0]->Fill(S.p);
                H1_probe_pT[0]->Fill(S.pT);
                H1_probe_theta[0]->Fill(S.theta);
                H1_probe_phi[0]->Fill(S.phi);
                H1_Q2[0]->Fill(K.Q2);
                H1_W2[0]->Fill(K.W2);
                H1_xB[0]->Fill(K.xB);
                H1_nu[0]->Fill(K.nu);
            }
        }
        nelectrons += Electrons.size();
        nphotons += Photons.size();
        H1_nelectrons[0]->Fill(Electrons.size());
        H1_nphotons[0]->Fill(Photons.size());
        H1_mon_nprobes->Fill(0.0, (int) Electrons.size() + Photons.size());
        //////////////////////
        //  AHDC::kftrack
        //////////////////////
        std::vector<State> Tracks;
        H1_ntracks[0]->Fill(trackBank.getRows()); 
        H1_mon_ntracks->Fill(0.0, (int) trackBank.getRows());
        for (int i = 0; i < trackBank.getRows(); i++) {
            State T(trackBank.getFloat("px",i), trackBank.getFloat("py",i), trackBank.getFloat("pz",i), 
                    trackBank.getFloat("x",i), trackBank.getFloat("y",i), trackBank.getFloat("z",i));
            T.vz = T.vz*0.1; // convert mm to cm
            T.p  = T.p*0.001; // convert MeV to GeV
            T.pT = T.pT*0.001; // convert MeV to Gev
            T.dEdx = track0Bank.getFloat("dEdx", i);
            T.adc = trackBank.getInt("sum_adc", i);
            T.trackid = trackBank.getInt("trackid", i);
            Tracks.push_back(T);
            // Histograms -- Level 0, no cuts on tracks
            H1_track_vz[0]->Fill(T.vz); 
            H1_track_p[0]->Fill(T.p); 
            H1_track_pT[0]->Fill(T.pT);
            H1_track_theta[0]->Fill(T.theta);
            H1_track_phi[0]->Fill(T.phi);
            H1_track_nhits[0]->Fill(trackBank.getInt("n_hits", i)); 
            H1_track_adc[0]->Fill(trackBank.getInt("sum_adc", i));
            H1_track_path[0]->Fill(trackBank.getFloat("path", i));   
            H1_track_dEdx[0]->Fill(track0Bank.getFloat("dEdx", i));    
            H1_track_residuals[0]->Fill(trackBank.getFloat("sum_residuals", i));
            H1_track_residuals_per_nhits[0]->Fill(trackBank.getFloat("sum_residuals", i)/trackBank.getInt("n_hits", i));
            H1_track_chi2[0]->Fill(trackBank.getFloat("chi2", i));
            H1_track_chi2_per_nhits[0]->Fill(trackBank.getFloat("chi2", i)/trackBank.getInt("n_hits", i));
            H1_track_p_drift[0]->Fill(0.001*trackBank.getFloat("p_drift", i)); // convert MeV to GeV
        }
        // Correlations -- level 0, no cuts
        // loop over all configurations
        for (State T : Tracks) {
            for (State S : Electrons) {
                H2_vze_vz[0]->Fill(S.vz, T.vz); 
                H1_delta_vz[0]->Fill(S.vz - T.vz); 
                H1_delta_phi[0]->Fill(S.phi - T.phi); // convert rad in deg 
                H2_p_dEdx[0]->Fill(S.pT, T.dEdx);
                H2_p_adc[0]->Fill(S.pT, T.adc);
                H2_pTe_pT[0]->Fill(S.pT, T.pT);
            }
            for (State S : Photons) {
                H2_vze_vz[0]->Fill(S.vz, T.vz); 
                H1_delta_vz[0]->Fill(S.vz - T.vz); 
                H1_delta_phi[0]->Fill(S.phi - T.phi); // convert rad in deg 
                H2_p_dEdx[0]->Fill(S.pT, T.dEdx);
                H2_p_adc[0]->Fill(S.pT, T.adc);
                H2_pTe_pT[0]->Fill(S.pT, T.pT);
            }
        }

        /*******************************************
         *  Define a priority order for the probes (cut level 1)
         * *****************************************/
        if (Electrons.size() > 0) {
            // select all electrons if there are; ignore photons
            H1_nelectrons[1]->Fill(Electrons.size());
            H1_mon_nprobes->Fill(1, Electrons.size());
            for (int i = 0; i < (int) Electrons.size(); i++) {
                State S = Electrons[i];
                Physics K = ElectronKinematics[i];
                H1_probe_vz[1]->Fill(S.vz);
                H1_probe_p[1]->Fill(S.p);
                H1_probe_pT[1]->Fill(S.pT);
                H1_probe_theta[1]->Fill(S.theta);
                H1_probe_phi[1]->Fill(S.phi);
                H1_Q2[1]->Fill(K.Q2);
                H1_W2[1]->Fill(K.W2);
                H1_xB[1]->Fill(K.xB);
                H1_nu[1]->Fill(K.nu);
            }
            Photons.clear();
            PhotonKinematics.clear();
        }
        else if (Photons.size() > 0) {
            std::vector<State> NewPhotons;
            std::vector<Physics> NewPhotonKinematics;
            for (int i = 0; i < (int) Photons.size(); i++) {
                State S = Photons[i];
                Physics K = PhotonKinematics[i];
                if (S.theta <= 5) { // keep only photons in the FT
                    H1_probe_vz[1]->Fill(S.vz);
                    H1_probe_p[1]->Fill(S.p);
                    H1_probe_pT[1]->Fill(S.pT);
                    H1_probe_theta[1]->Fill(S.theta);
                    H1_probe_phi[1]->Fill(S.phi);
                    H1_Q2[1]->Fill(K.Q2);
                    H1_W2[1]->Fill(K.W2);
                    H1_xB[1]->Fill(K.xB);
                    H1_nu[1]->Fill(K.nu);
                    NewPhotons.push_back(S);
                    NewPhotonKinematics.push_back(K);
                }
            }
            H1_nphotons[1]->Fill(NewPhotons.size());
            H1_mon_nprobes->Fill(1, NewPhotons.size());
            // keep only photons in the FT
            Photons = NewPhotons;
            PhotonKinematics = NewPhotonKinematics;
        }
        else {
            continue; // ignore this event
        }
        // if the event is not ignored, let check the tracks
        Tracks.clear(); // restart
        H1_ntracks[1]->Fill(trackBank.getRows()); 
        H1_mon_ntracks->Fill(1.0, (int) trackBank.getRows());
        for (int i = 0; i < trackBank.getRows(); i++) {
            State T(trackBank.getFloat("px",i), trackBank.getFloat("py",i), trackBank.getFloat("pz",i), 
                    trackBank.getFloat("x",i), trackBank.getFloat("y",i), trackBank.getFloat("z",i));
            T.vz = T.vz*0.1; // convert mm to cm
            T.p  = T.p*0.001; // convert MeV to GeV
            T.pT = T.pT*0.001; // convert MeV to Gev
            T.dEdx = track0Bank.getFloat("dEdx", i);
            T.adc = trackBank.getInt("sum_adc", i);
            T.trackid = trackBank.getInt("trackid", i);
            Tracks.push_back(T);
            // Histograms -- Level 1
            H1_track_vz[1]->Fill(T.vz);
            H1_track_p[1]->Fill(T.p); 
            H1_track_pT[1]->Fill(T.pT);
            H1_track_theta[1]->Fill(T.theta);
            H1_track_phi[1]->Fill(T.phi);
            H1_track_nhits[1]->Fill(trackBank.getInt("n_hits", i)); 
            H1_track_adc[1]->Fill(trackBank.getInt("sum_adc", i));
            H1_track_path[1]->Fill(trackBank.getFloat("path", i));   
            H1_track_dEdx[1]->Fill(track0Bank.getFloat("dEdx", i));    
            H1_track_residuals[1]->Fill(trackBank.getFloat("sum_residuals", i));
            H1_track_residuals_per_nhits[1]->Fill(trackBank.getFloat("sum_residuals", i)/trackBank.getInt("n_hits", i));
            H1_track_chi2[1]->Fill(trackBank.getFloat("chi2", i));
            H1_track_chi2_per_nhits[1]->Fill(trackBank.getFloat("chi2", i)/trackBank.getInt("n_hits", i));
            H1_track_p_drift[1]->Fill(0.001*trackBank.getFloat("p_drift", i)); // convert MeV to GeV
        }
        // Now Combined Electrons and Photons
        std::vector<State> Probes;
        Probes.reserve(Electrons.size() + Photons.size());
        Probes.insert(Probes.end(), Electrons.begin(), Electrons.end());
        Probes.insert(Probes.end(), Photons.begin(), Photons.end());
        std::vector<Physics> ProbeKinematics;
        ProbeKinematics.reserve(ElectronKinematics.size() + PhotonKinematics.size());
        ProbeKinematics.insert(ProbeKinematics.end(), ElectronKinematics.begin(), ElectronKinematics.end());
        ProbeKinematics.insert(ProbeKinematics.end(), PhotonKinematics.begin(), PhotonKinematics.end());
        // let check the correlations
        for (State T : Tracks) {
            for (State S : Probes) {
                H2_vze_vz[1]->Fill(S.vz, T.vz); 
                H1_delta_vz[1]->Fill(S.vz - T.vz); 
                H1_delta_phi[1]->Fill(S.phi - T.phi);  
                H2_p_dEdx[1]->Fill(S.pT, T.dEdx);
                H2_p_adc[1]->Fill(S.pT, T.adc);
                H2_pTe_pT[1]->Fill(S.pT, T.pT);
            }
        }
        /*******************************************
         *  cut on W2 (cut level 2)
         * *****************************************/
        int count_electrons2 = 0;
        int count_photons2 = 0;
        int count_tracks2 = 0;
        int count_electrons3 = 0; // level 3
        int count_photons3 = 0;
        int count_tracks3 = 0;
        std::vector<ElasticsOutput> Elastics; // list of couples (e/gamma, track) for elastics
        for (int i = 0; i < (int) Probes.size(); i++) {
            State S = Probes[i];
            Physics K = ProbeKinematics[i];
            if ((K.W2 < W2_min) || (K.W2 > W2_max)) continue;
            if (S.pid == 11) count_electrons2++;
            if (S.pid == 22) count_photons2++;
            for (State T : Tracks) {
                count_tracks2++;
                // corr
                H2_vze_vz[2]->Fill(S.vz, T.vz); 
                H1_delta_vz[2]->Fill(S.vz - T.vz); 
                H1_delta_phi[2]->Fill(S.phi - T.phi);  
                H2_p_dEdx[2]->Fill(S.pT, T.dEdx);
                H2_p_adc[2]->Fill(S.pT, T.adc);
                H2_pTe_pT[2]->Fill(S.pT, T.pT);
                // probe 
                H1_probe_vz[2]->Fill(S.vz);
                H1_probe_p[2]->Fill(S.p);
                H1_probe_pT[2]->Fill(S.pT);
                H1_probe_theta[2]->Fill(S.theta);
                H1_probe_phi[2]->Fill(S.phi);
                H1_Q2[2]->Fill(K.Q2);
                H1_W2[2]->Fill(K.W2);
                H1_xB[2]->Fill(K.xB);
                H1_nu[2]->Fill(K.nu);
                // track
                H1_track_vz[2]->Fill(T.vz);
                H1_track_p[2]->Fill(T.p);
                H1_track_pT[2]->Fill(T.pT);
                H1_track_theta[2]->Fill(T.theta);
                H1_track_phi[2]->Fill(T.phi);
                H1_track_nhits[2]->Fill(trackBank.getInt("n_hits", i)); 
                H1_track_adc[2]->Fill(trackBank.getInt("sum_adc", i));
                H1_track_path[2]->Fill(trackBank.getFloat("path", i));   
                H1_track_dEdx[2]->Fill(track0Bank.getFloat("dEdx", i));    
                H1_track_residuals[2]->Fill(trackBank.getFloat("sum_residuals", i));
                H1_track_residuals_per_nhits[2]->Fill(trackBank.getFloat("sum_residuals", i)/trackBank.getInt("n_hits", i));
                H1_track_chi2[2]->Fill(trackBank.getFloat("chi2", i));
                H1_track_chi2_per_nhits[2]->Fill(trackBank.getFloat("chi2", i)/trackBank.getInt("n_hits", i));
                H1_track_p_drift[2]->Fill(0.001*trackBank.getFloat("p_drift", i)); // convert MeV to GeV
                /*******************************************
                 *  delta_phi (cut level 3)
                 * *****************************************/
                double delta_phi = S.phi - T.phi;
                double stdev = 20;
                if ((fabs(delta_phi+198) < stdev) || (fabs(delta_phi-162) < stdev)) {
                    if (S.pid == 11) count_electrons3++;
                    if (S.pid == 22) count_photons3++;
                    count_tracks3++;
                    // corr
                    H2_vze_vz[3]->Fill(S.vz, T.vz); 
                    H1_delta_vz[3]->Fill(S.vz - T.vz); 
                    H1_delta_phi[3]->Fill(S.phi - T.phi);  
                    H2_p_dEdx[3]->Fill(S.pT, T.dEdx);
                    H2_p_adc[3]->Fill(S.pT, T.adc);
                    H2_pTe_pT[3]->Fill(S.pT, T.pT);
                    // probe 
                    H1_probe_vz[3]->Fill(S.vz);
                    H1_probe_p[3]->Fill(S.p);
                    H1_probe_pT[3]->Fill(S.pT);
                    H1_probe_theta[3]->Fill(S.theta);
                    H1_probe_phi[3]->Fill(S.phi);
                    H1_Q2[3]->Fill(K.Q2);
                    H1_W2[3]->Fill(K.W2);
                    H1_xB[3]->Fill(K.xB);
                    H1_nu[3]->Fill(K.nu);
                    // track
                    H1_track_vz[3]->Fill(T.vz);
                    H1_track_p[3]->Fill(T.p); 
                    H1_track_pT[3]->Fill(T.pT);
                    H1_track_theta[3]->Fill(T.theta);
                    H1_track_phi[3]->Fill(T.phi);
                    H1_track_nhits[3]->Fill(trackBank.getInt("n_hits", i)); 
                    H1_track_adc[3]->Fill(trackBank.getInt("sum_adc", i));
                    H1_track_path[3]->Fill(trackBank.getFloat("path", i));   
                    H1_track_dEdx[3]->Fill(track0Bank.getFloat("dEdx", i));    
                    H1_track_residuals[3]->Fill(trackBank.getFloat("sum_residuals", i));
                    H1_track_residuals_per_nhits[3]->Fill(trackBank.getFloat("sum_residuals", i)/trackBank.getInt("n_hits", i));
                    H1_track_chi2[3]->Fill(trackBank.getFloat("chi2", i));
                    H1_track_chi2_per_nhits[3]->Fill(trackBank.getFloat("chi2", i)/trackBank.getInt("n_hits", i));
                    H1_track_p_drift[3]->Fill(0.001*trackBank.getFloat("p_drift", i)); // convert MeV to GeV
                    /***********************************************
                     * for simulation calibration
                     * ********************************************/
                    Elastics.push_back({S,T}); 
                }
            }
        }
        // level 2
        H1_nelectrons[2]->Fill(count_electrons2);
        H1_nphotons[2]->Fill(count_photons2);
        H1_ntracks[2]->Fill(count_tracks2);
        H1_mon_nprobes->Fill(2.0, count_electrons2 + count_photons2);
        H1_mon_ntracks->Fill(2.0, count_tracks2);
        // level 3
        H1_nelectrons[3]->Fill(count_electrons3);
        H1_nphotons[3]->Fill(count_photons3);
        H1_ntracks[3]->Fill(count_tracks3);
        H1_mon_nprobes->Fill(3.0, count_electrons3 + count_photons3);
        H1_mon_ntracks->Fill(3.0, count_tracks3);
        /***********************************************
         * for simulation calibration
         * ********************************************/
        std::vector<int> AdcId;
        std::vector<int> TrackId;
        H1_nelastics->Fill(Elastics.size());
        for (ElasticsOutput e : Elastics) {
            TrackId.push_back(e.track.trackid); 
            for (int i = 0; i < hitBank.getRows(); i++) {
                int trackid = hitBank.getInt("trackid",i);
                if (trackid == e.track.trackid) {
                    // hits
                    int hit_id = hitBank.getInt("id",i); // they are the column index in AHDC::adc bank
                    double sector    = adcBank.getInt("sector", hit_id);
                    double layer     = adcBank.getInt("layer", hit_id);
                    double component = adcBank.getInt("component", hit_id);
                    ahdcT0 obj  = ahdcConstants.get_t0(sector, layer, component);
                    double time    = adcBank.getFloat("leadingEdgeTime", hit_id);
                    double timeMax = adcBank.getFloat("time", hit_id);
                    double tot     = adcBank.getFloat("timeOverThreshold", hit_id);
                    int adc        = adcBank.getInt("ADC", hit_id);
                    // hit cut on time
                    H1_dt00->Fill(time-obj.t0);
                    if ((time - obj.t0 < 50) || (time - obj.t0 > 250)) continue;
                    H1_tot0->Fill(tot);
                    if ((tot < 400) || (tot > 600)) continue;
                    //int wfType = adcBank.getInt("wfType", hit_id);
                    //H1_wfType->Fill(wfType);
                    H1_t0->Fill(obj.t0);
                    H1_dt0->Fill(time-obj.t0);
                    H1_leadingEdgeTime->Fill(time); 
                    H1_timeMax->Fill(timeMax);
                    H1_deltaTime->Fill(timeMax-time); 
                    H1_timeOverThreshold->Fill(tot);
                    H1_amplitude->Fill(adc); 
                    H2_times->Fill(timeMax, time);
                    H2_amp_tot->Fill(adc, tot);
                    H2_deltaTime_adc->Fill(timeMax-time, adc); 
                    // track and probe
                    H1_elastic_track_p->Fill(e.track.p);
                    H1_elastic_track_pT->Fill(e.track.pT);
                    H1_elastic_track_theta->Fill(e.track.theta);
                    H1_elastic_track_phi->Fill(e.track.phi);
                    H1_elastic_probe_p->Fill(e.probe.p);
                    H1_elastic_probe_pT->Fill(e.probe.pT);
                    H1_elastic_probe_theta->Fill(e.probe.theta);
                    H1_elastic_probe_phi->Fill(e.probe.phi);
                    H2_elastic_pT_adc->Fill(e.probe.pT, e.track.adc);
                    // theoretical values of the track given the electron kinematics
                    double theta = M_PI*e.probe.theta/180;
                    double p = 2*Ee*sin(theta/2);
                    double pT = Ee*sin(theta);
                    double pz = Ee*(1-cos(theta));
                    double theta_track = acos(pz/p)*180/M_PI;
                    double phi_track = (e.probe.phi > 180) ? (e.probe.phi - 180) : (e.probe.phi + 180);
                    H1_elastic_theory_p->Fill(p);
                    H1_elastic_theory_pT->Fill(pT);
                    H1_elastic_theory_theta->Fill(theta_track);
                    H1_elastic_theory_phi->Fill(phi_track);
                    if ((e.probe.pT > 0.230) && (e.probe.pT < 0.260) && (e.track.adc > 1000) && (e.track.adc < 3500)) {
                        AdcId.push_back(hit_id); // they correspond to the column index in AHDC::adc
                        H1_selection_theory_p->Fill(p);
                        H1_selection_theory_pT->Fill(pT);
                        H1_selection_theory_theta->Fill(theta_track);
                        H1_selection_theory_phi->Fill(phi_track);
                        H1_selection_theory_adc->Fill(e.track.adc);
                    }
                }
            } // end loop over elastics
        } // end loop over elastics hits
        // save those event in a hipo file
        if (AdcId.size()*TrackId.size() > 0) { 
            hipo::bank outAdcBank(schemaAdc, (int) AdcId.size());
            hipo::bank outWfBank(schemaWf, (int) AdcId.size());
            hipo::bank outTrackBank(schemaTrack, (int) TrackId.size());
            for (int i = 0; i < outAdcBank.getRows(); i++) {
                // Attention: the hit id numerotation starts at 1, the getter indexation starts at 0
                int index = AdcId[i] - 1;
                // AHDC::adc
                outAdcBank.putInt("sector", i, adcBank.getInt("sector", index));
                outAdcBank.putInt("layer", i, adcBank.getInt("layer", index));
                outAdcBank.putInt("component", i, adcBank.getInt("component", index));
                outAdcBank.putInt("order", i, adcBank.getInt("order", index));
                outAdcBank.putInt("ADC", i, adcBank.getInt("ADC", index));
                outAdcBank.putFloat("time", i, adcBank.getFloat("time", index));
                outAdcBank.putFloat("ped", i, adcBank.getInt("ped", index));
                outAdcBank.putInt("windex", i, adcBank.getInt("windex", index));
                outAdcBank.putInt("integral", i, adcBank.getInt("integral", index));
                outAdcBank.putFloat("leadingEdgeTime", i, adcBank.getFloat("leadingEdgeTime", index));
                outAdcBank.putFloat("timeOverThreshold", i, adcBank.getFloat("timeOverThreshold", index));
                outAdcBank.putFloat("constantFractionTime", i, adcBank.getFloat("constantFractionTime", index));
                //outAdcBank.put("wfType", i, adcBank.getInt("ped", index));
                outAdcBank.putInt("wfType", i, 0);
                // AHDC::wf
                outWfBank.putInt("sector", i, wfBank.getInt("sector", index));
                outWfBank.putInt("layer", i, wfBank.getInt("layer", index));
                outWfBank.putInt("component", i, wfBank.getInt("component", index));
                outWfBank.putLong("timestamp", i, wfBank.getLong("timestamp", index));
                outWfBank.putInt("time", i, wfBank.getInt("time", index));
                outWfBank.putInt("s1", i, wfBank.getInt("s1", index));
                outWfBank.putInt("s2", i, wfBank.getInt("s2", index));
                outWfBank.putInt("s3", i, wfBank.getInt("s3", index));
                outWfBank.putInt("s4", i, wfBank.getInt("s4", index));
                outWfBank.putInt("s5", i, wfBank.getInt("s5", index));
                outWfBank.putInt("s6", i, wfBank.getInt("s6", index));
                outWfBank.putInt("s7", i, wfBank.getInt("s7", index));
                outWfBank.putInt("s8", i, wfBank.getInt("s8", index));
                outWfBank.putInt("s9", i, wfBank.getInt("s9", index));
                outWfBank.putInt("s10", i, wfBank.getInt("s10", index));
                outWfBank.putInt("s11", i, wfBank.getInt("s11", index));
                outWfBank.putInt("s12", i, wfBank.getInt("s12", index));
                outWfBank.putInt("s13", i, wfBank.getInt("s13", index));
                outWfBank.putInt("s14", i, wfBank.getInt("s14", index));
                outWfBank.putInt("s15", i, wfBank.getInt("s15", index));
                outWfBank.putInt("s16", i, wfBank.getInt("s16", index));
                outWfBank.putInt("s17", i, wfBank.getInt("s17", index));
                outWfBank.putInt("s18", i, wfBank.getInt("s18", index));
                outWfBank.putInt("s19", i, wfBank.getInt("s19", index));
                outWfBank.putInt("s20", i, wfBank.getInt("s20", index));
                outWfBank.putInt("s21", i, wfBank.getInt("s21", index));
                outWfBank.putInt("s22", i, wfBank.getInt("s22", index));
                outWfBank.putInt("s23", i, wfBank.getInt("s23", index));
                outWfBank.putInt("s24", i, wfBank.getInt("s24", index));
                outWfBank.putInt("s25", i, wfBank.getInt("s25", index));
                outWfBank.putInt("s26", i, wfBank.getInt("s26", index));
                outWfBank.putInt("s27", i, wfBank.getInt("s27", index));
                outWfBank.putInt("s28", i, wfBank.getInt("s28", index));
                outWfBank.putInt("s29", i, wfBank.getInt("s29", index));
                outWfBank.putInt("s30", i, wfBank.getInt("s30", index));
            } 
            for (int i = 0; i < outTrackBank.getRows(); i++) {
                // AHDC::track
                // Attention: the track id numerotation starts at 1, the getter indexation starts at 0
                int index = TrackId[i] - 1;
                //printf("trackid: %d/%d ", TrackId[i], index);
                outTrackBank.putInt("trackid", i, trackBank.getInt("trackid", index));
                outTrackBank.putFloat("x", i, trackBank.getFloat("x", index));
                outTrackBank.putFloat("y", i, trackBank.getFloat("y", index));
                outTrackBank.putFloat("z", i, trackBank.getFloat("z", index));
                outTrackBank.putFloat("px", i, trackBank.getFloat("px", index));
                outTrackBank.putFloat("py", i, trackBank.getFloat("py", index));
                //printf("pz: %lf ", trackBank.getFloat("pz", index));
                outTrackBank.putFloat("pz", i, trackBank.getFloat("pz", index));
                outTrackBank.putInt("n_hits", i, trackBank.getInt("n_hits", index));
                //printf("sum_adc: %d ", trackBank.getInt("sum_adc", index));
                outTrackBank.putInt("sum_adc", i, trackBank.getInt("sum_adc", index));
                outTrackBank.putFloat("path", i, trackBank.getFloat("path", index));
                //printf("path: %lf ", trackBank.getFloat("path", index));
                outTrackBank.putFloat("dEdx", i, trackBank.getFloat("dEdx", index));
                outTrackBank.putFloat("p_drift", i, trackBank.getFloat("p_drift", index));
                outTrackBank.putFloat("chi2", i, trackBank.getFloat("chi2", index));
                outTrackBank.putFloat("sum_residuals", i, trackBank.getFloat("sum_residuals", index));
            } // end write selection from elastics in a hipo file
            outEvent.reset();
            outEvent.addStructure(outAdcBank);
            outEvent.addStructure(outWfBank);
            outEvent.addStructure(outTrackBank);
            writer.addEvent(outEvent);
        }
    } // and loop over events
    writer.close();
    printf("Create file : %s\n", outputFile);

    printf("nevents    : %ld \n", nevents);
    printf("nelectrons : %ld \n", nelectrons);
    printf("nphotons   : %ld \n", nphotons);
    printf("ntracks    : %ld \n", ntracks);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("* time elapsed : %lf seconds\n", elapsed.count());

/************************************************************
 * output
 * **********************************************************/

    const char * output = "./output/corrected_elastics.root";
    TFile *f = new TFile(output, "RECREATE");
    TDirectory *probes_dir   = f->mkdir("electrons_photons");
    TDirectory *physics_dir  = f->mkdir("physics");
    TDirectory *tracks_dir  = f->mkdir("tracks");
    TDirectory *corr_dir     = f->mkdir("correlations");
    TDirectory *calibration     = f->mkdir("calibration");
    
    //////////////////////
    //  probes
    //////////////////////
    f->cd();
    H1_mon_nprobes->Write("nprobes_evolution");
    H1_mon_ntracks->Write("ntracks_evolution");
    // Level 0 -- no cuts
    TDirectory* probes_nocuts_dir = probes_dir->mkdir("nocuts");
    probes_nocuts_dir->cd();
    H1_probe_p[0]->Write("p");
    H1_probe_pT[0]->Write("pT");
    H1_probe_vz[0]->Write("vz");
    H1_probe_theta[0]->Write("theta");
    H1_probe_phi[0]->Write("phi");
    H1_nelectrons[0]->Write("nelectrons");
    H1_nphotons[0]->Write("nphotons");
    // Level 1 -- priority on the selection of the probe
    TDirectory* probes_priority_dir = probes_dir->mkdir("priority");
    probes_priority_dir->cd();
    H1_probe_p[1]->Write("p");
    H1_probe_pT[1]->Write("pT");
    H1_probe_vz[1]->Write("vz");
    H1_probe_theta[1]->Write("theta");
    H1_probe_phi[1]->Write("phi");
    H1_nelectrons[1]->Write("nelectrons");
    H1_nphotons[1]->Write("nphotons");
    // Level 2 -- cut_on_W2
    TDirectory* probes_delta_phi_dir = probes_dir->mkdir("cut_on_W2");
    probes_delta_phi_dir->cd();
    H1_probe_p[2]->Write("p");
    H1_probe_pT[2]->Write("pT");
    H1_probe_vz[2]->Write("vz");
    H1_probe_theta[2]->Write("theta");
    H1_probe_phi[2]->Write("phi");
    H1_nelectrons[2]->Write("nelectrons");
    H1_nphotons[2]->Write("nphotons");
    // Level 3 -- cut_on_delta_phi
    TDirectory* probes_w2_cut_dir = probes_dir->mkdir("cut_on_delta_phi");
    probes_w2_cut_dir->cd();
    H1_probe_p[3]->Write("p");
    H1_probe_pT[3]->Write("pT");
    H1_probe_vz[3]->Write("vz");
    H1_probe_theta[3]->Write("theta");
    H1_probe_phi[3]->Write("phi");
    H1_nelectrons[3]->Write("nelectrons");
    H1_nphotons[3]->Write("nphotons");
    //////////////////////
    //  physics
    //////////////////////
    // Level 0 -- no cuts
    TDirectory* physics_nocuts_dir = physics_dir->mkdir("nocuts");
    physics_nocuts_dir->cd();
    H1_W2[0]->Write("W2");
    H1_Q2[0]->Write("Q2");
    H1_xB[0]->Write("xB");
    H1_nu[0]->Write("nu");
    // Level 1 -- priority on the selection of the probe
    TDirectory* physics_priority_dir = physics_dir->mkdir("priority");
    physics_priority_dir->cd();
    H1_W2[1]->Write("W2");
    H1_Q2[1]->Write("Q2");
    H1_xB[1]->Write("xB");
    H1_nu[1]->Write("nu");
    // Level 2 -- cut_on_W2
    TDirectory* physics_delta_phi_dir = physics_dir->mkdir("cut_on_W2");
    physics_delta_phi_dir->cd();
    H1_W2[2]->Write("W2");
    H1_Q2[2]->Write("Q2");
    H1_xB[2]->Write("xB");
    H1_nu[2]->Write("nu");
    // Level 3 -- cut_on_delta_phi
    TDirectory* physics_w2_cut_dir = physics_dir->mkdir("cut_on_delta_phi");
    physics_w2_cut_dir->cd();
    H1_W2[3]->Write("W2");
    H1_Q2[3]->Write("Q2");
    H1_xB[3]->Write("xB");
    H1_nu[3]->Write("nu");
    //////////////////////
    //  tracks
    //////////////////////
    // Level 0 -- no cuts
    TDirectory* tracks_nocuts_dir = tracks_dir->mkdir("nocuts");
    tracks_nocuts_dir->cd();
    H1_ntracks[0]->Write("ntracks");
    H1_track_vz[0]->Write("vz");
    H1_track_p[0]->Write("p");
    H1_track_p_drift[0]->Write("p_drift");
    H1_track_pT[0]->Write("pT");
    H1_track_theta[0]->Write("theta");
    H1_track_phi[0]->Write("phi");
    H1_track_nhits[0]->Write("nhits");
    H1_track_adc[0]->Write("adc");
    H1_track_path[0]->Write("path");
    H1_track_dEdx[0]->Write("dEdx");
    H1_track_residuals[0]->Write("residuals");
    H1_track_residuals_per_nhits[0]->Write("residuals_per_hit");
    H1_track_chi2[0]->Write("chi2");
    H1_track_chi2_per_nhits[0]->Write("chi2_per_hit");
    // Level 1 -- probe priority
    TDirectory* tracks_probe_dir = tracks_dir->mkdir("probe_priority");
    tracks_probe_dir->cd();
    H1_ntracks[1]->Write("ntracks");
    H1_track_vz[1]->Write("vz");
    H1_track_p[1]->Write("p");
    H1_track_p_drift[1]->Write("p_drift");
    H1_track_pT[1]->Write("pT");
    H1_track_theta[1]->Write("theta");
    H1_track_phi[1]->Write("phi");
    H1_track_nhits[1]->Write("nhits");
    H1_track_adc[1]->Write("adc");
    H1_track_path[1]->Write("path");
    H1_track_dEdx[1]->Write("dEdx");
    H1_track_residuals[1]->Write("residuals");
    H1_track_residuals_per_nhits[1]->Write("residuals_per_hit");
    H1_track_chi2[1]->Write("chi2");
    H1_track_chi2_per_nhits[1]->Write("chi2_per_hit");
    // Level 2 -- cut_on_W2
    TDirectory* tracks_delta_phi_dir = tracks_dir->mkdir("cut_on_W2");
    tracks_delta_phi_dir->cd();
    H1_ntracks[2]->Write("ntracks");
    H1_track_vz[2]->Write("vz");
    H1_track_p[2]->Write("p");
    H1_track_p_drift[2]->Write("p_drift");
    H1_track_pT[2]->Write("pT");
    H1_track_theta[2]->Write("theta");
    H1_track_phi[2]->Write("phi");
    H1_track_nhits[2]->Write("nhits");
    H1_track_adc[2]->Write("adc");
    H1_track_path[2]->Write("path");
    H1_track_dEdx[2]->Write("dEdx");
    H1_track_residuals[2]->Write("residuals");
    H1_track_residuals_per_nhits[2]->Write("residuals_per_hit");
    H1_track_chi2[2]->Write("chi2");
    H1_track_chi2_per_nhits[2]->Write("chi2_per_hit");
    // Level 3 -- cut_on_delta_phi
    TDirectory* tracks_w2_cut_dir = tracks_dir->mkdir("cut_on_delta_phi");
    tracks_w2_cut_dir->cd();
    H1_ntracks[3]->Write("ntracks");
    H1_track_vz[3]->Write("vz");
    H1_track_p[3]->Write("p");
    H1_track_p_drift[3]->Write("p_drift");
    H1_track_pT[3]->Write("pT");
    H1_track_theta[3]->Write("theta");
    H1_track_phi[3]->Write("phi");
    H1_track_nhits[3]->Write("nhits");
    H1_track_adc[3]->Write("adc");
    H1_track_path[3]->Write("path");
    H1_track_dEdx[3]->Write("dEdx");
    H1_track_residuals[3]->Write("residuals");
    H1_track_residuals_per_nhits[3]->Write("residuals_per_hit");
    H1_track_chi2[3]->Write("chi2");
    H1_track_chi2_per_nhits[3]->Write("chi2_per_hit");
    //////////////////////
    //  correlation
    //////////////////////
    for (TH2D* ptr : H2_vze_vz) {
        ptr->GetXaxis()->SetTitle("vz probe (cm)");
        ptr->GetYaxis()->SetTitle("vz track (cm)");
    }
    for (TH2D* ptr : H2_pTe_pT) {
        ptr->GetXaxis()->SetTitle("pT probe (GeV)");
        ptr->GetYaxis()->SetTitle("pT track (GeV)");
    }
    for (TH2D* ptr : H2_p_dEdx) {
        ptr->GetXaxis()->SetTitle("pT probe (GeV)");
        ptr->GetYaxis()->SetTitle("dEdx track (adc/mm)");
    }
    for (TH2D* ptr : H2_p_adc) {
        ptr->GetXaxis()->SetTitle("pT probe (GeV)");
        ptr->GetYaxis()->SetTitle("adc track (adc)");
    }
    // Level 0 -- no cuts
    TDirectory* corr_nocuts_dir = corr_dir->mkdir("nocuts");
    corr_nocuts_dir->cd();
    H2_vze_vz[0]->Write("corr_vz"); 
    H1_delta_vz[0]->Write("delta_vz"); 
    H1_delta_phi[0]->Write("delat_phi");  
    H2_p_dEdx[0]->Write("corr_pT_dEdx");
    H2_p_adc[0]->Write("corr_pT_adc");
    H2_pTe_pT[0]->Write("corr_pT");
    // Level 1 -- probe priority
    TDirectory* corr_probe_priority_dir = corr_dir->mkdir("probe_priority");
    corr_probe_priority_dir->cd();
    H2_vze_vz[1]->Write("corr_vz"); 
    H1_delta_vz[1]->Write("delta_vz"); 
    H1_delta_phi[1]->Write("delat_phi");  
    H2_p_dEdx[1]->Write("corr_pT_dEdx");
    H2_p_adc[1]->Write("corr_pT_adc");
    H2_pTe_pT[1]->Write("corr_pT");
    // Level 2 -- cu_on_W2
    TDirectory* corr_delta_phi_dir = corr_dir->mkdir("cut_on_W2");
    corr_delta_phi_dir->cd();
    H2_vze_vz[2]->Write("corr_vz"); 
    H1_delta_vz[2]->Write("delta_vz"); 
    H1_delta_phi[2]->Write("delat_phi");  
    H2_p_dEdx[2]->Write("corr_pT_dEdx");
    H2_p_adc[2]->Write("corr_pT_adc");
    H2_pTe_pT[2]->Write("corr_pT");
    // Level 3 -- cut_on_delta_phi
    TDirectory* corr_w2_cut_dir = corr_dir->mkdir("cut_on_delta_phi");
    corr_w2_cut_dir->cd();
    H2_vze_vz[3]->Write("corr_vz"); 
    H1_delta_vz[3]->Write("delta_vz"); 
    H1_delta_phi[3]->Write("delat_phi");  
    H2_p_dEdx[3]->Write("corr_pT_dEdx");
    H2_p_adc[3]->Write("corr_pT_adc");
    H2_pTe_pT[3]->Write("corr_pT");
    //////////////////////
    //  calibration
    //////////////////////
    calibration->cd(); 
    H1_dt00->Write("dt00"); 
    H1_dt0->Write("dt0"); 
    H1_tot0->Write("tot0"); 
    H1_timeOverThreshold->Write("tot");
    H1_leadingEdgeTime->Write("time"); 
    H1_timeMax->Write("timeMax");
    H1_t0->Write("t0"); 
    H1_deltaTime->Write("dt"); 
    H1_wfType->Write("wfType");
    H1_amplitude->Write("amplitude"); 
    H2_times->Write("corr_times");
    H2_amp_tot->Write("corr_amp_tot");
    H2_deltaTime_adc->Write("corr_dt_amp");
    // >>>>>>>>>   Process H2_deltaTime_adc
    TCanvas* canvas1 = new TCanvas("c1_deltaTime_amplitude", "deltaTime vs amplitude; deltaTime (ns); amplitude (adc)");
    canvas1->SetLogz();
    H2_deltaTime_adc->Draw("colz");
    int nbins = H2_deltaTime_adc->GetXaxis()->GetNbins();
    printf("> Analyse correlation between deltaTime (ns) and amplitude (adc)\n");
    printf(">   nbins x axis : %d\n", nbins);
    int proj_nbins = 0.05*nbins;
    int Npts = nbins/proj_nbins;
    TGraphErrors* gr = new TGraphErrors(Npts);
    for (int i = 0; i < Npts; i++) {
        TH1D* H1_projTime = H2_deltaTime_adc->ProjectionX(TString::Format("_px_%d", i).Data(), i*proj_nbins, (i+1)*proj_nbins);
        double mean = H1_projTime->GetMean();
        double stdDev = H1_projTime->GetStdDev();
        printf(">   mean : %lf, stdDev : %lf\n", mean, stdDev);
        double y_bin_inf = H2_deltaTime_adc->GetYaxis()->GetBinCenter(i*proj_nbins); // adc_inf for this collection of bins
        double y_bin_sup = H2_deltaTime_adc->GetYaxis()->GetBinCenter((i+1)*proj_nbins); // adc_sup this collection of bins
        gr->SetPoint(i, mean, 0.5*(y_bin_inf+y_bin_sup));
        gr->SetPointError(i, stdDev, 0);
    }
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(21);
    gr->SetLineWidth(2);
    gr->SetLineColor(kRed);
    gr->Draw("same pl");
    TLegend* legend = new TLegend(0.1,0.7,0.25,0.9);
    //legend->SetHeader("deltaTime projection","C");
    legend->AddEntry(gr, "deltaTime projection");
    legend->Draw();
    // <<<<<<<<<<<<<< end process H2_deltaTime_adc
    canvas1->Write("c1_deltaTime_adc");
    H1_elastic_track_p->Write("e.track.p");
    H1_elastic_track_pT->Write("e.track.pT");
    H1_elastic_track_theta->Write("e.track.theta");
    H1_elastic_track_phi->Write("e.track.phi");
    H1_elastic_probe_p->Write("e.probe.p");
    H1_elastic_probe_pT->Write("e.probe.pT");
    H1_elastic_probe_theta->Write("e.probe.theta");
    H1_elastic_probe_phi->Write("e.probe.phi");
    H2_elastic_pT_adc->Write("e.pT_vs_adc");
    H1_photon_p->Write("ft_photon_p");
    H1_photon_pT->Write("ft_photon_pT");
    H1_photon_theta->Write("ft_photon_theta");
    H1_photon_phi->Write("ft_photon_phi");
    H1_nelastics->Write("nelastics");
    H1_elastic_theory_p->Write("theory.track.p");
    H1_elastic_theory_pT->Write("theory.track.pT");
    H1_elastic_theory_theta->Write("theory.track.theta");
    H1_elastic_theory_phi->Write("theory.track.phi");
    TCanvas* canvas2 = new TCanvas("c2_elastics", "real track vs theory");
    canvas2->Divide(2,2);
    canvas2->cd(1);
    THStack* stack1 = new THStack("stack_p", "p #color[4]{reconstructed} vs #color[2]{theory}; p (GeV)");
    //TLegend* legend2 = new TLegend(0.1,0.7,0.35,0.9);
    //legend->AddEntry(H1_elastic_track_p, "reconstructed track");
    stack1->Add(H1_elastic_track_p);
    H1_elastic_theory_p->SetLineColor(kRed);
    stack1->Add(H1_elastic_theory_p);
    //legend->AddEntry(H1_elastic_theory_p, "theoretical track");
    //legend2->Draw();
    stack1->Draw("nostack");
    canvas2->cd(2);
    THStack* stack2 = new THStack("stack_pT", "pT #color[4]{reconstructed} vs #color[2]{theory}; pT (GeV)");
    stack2->Add(H1_elastic_track_pT);
    H1_elastic_theory_pT->SetLineColor(kRed);
    stack2->Add(H1_elastic_theory_pT);
    stack2->Draw("nostack");
    canvas2->cd(3);
    THStack* stack3 = new THStack("stack_theta", "theta #color[4]{reconstructed} vs #color[2]{theory}; theta (deg)");
    stack3->Add(H1_elastic_track_theta);
    H1_elastic_theory_theta->SetLineColor(kRed);
    stack3->Add(H1_elastic_theory_theta);
    stack3->Draw("nostack");
    canvas2->cd(4);
    THStack* stack4 = new THStack("stack_phi", "phi #color[4]{reconstructed} vs #color[2]{theory}; phi (deg)");
    stack4->Add(H1_elastic_track_phi);
    H1_elastic_theory_phi->SetLineColor(kRed);
    stack4->Add(H1_elastic_theory_phi);
    stack4->Draw("nostack");
    canvas2->Write("do_real_track_theory?");
    // select proton
    H1_selection_theory_p->Write("selection.theory.track.p");
    H1_selection_theory_pT->Write("selection.theory.track.pT");
    H1_selection_theory_theta->Write("selection.theory.track.theta");
    H1_selection_theory_phi->Write("selection_.theory.track.phi");
    H1_selection_theory_adc->Write("selection_.theory.adc.phi");

    //////////////////////////////
    // Close file
    // ///////////////////////////
    f->Close();
    printf("File created : %s\n", output);

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

State::State(double _px, double _py, double _pz, double _vx, double _vy, double _vz) 
    : px(_px), py(_py), pz(_pz), vx(_vx), vy(_vy), vz(_vz) {
    futils::cart2polar(px, py, pz, p, theta, phi);
    pT = sqrt(px*px + py*py);
    // convert theta and phi in deg
    theta = theta*180/M_PI;
    phi = phi*180/M_PI;
}

Physics::Physics(State Probe, double _Mt, double _Ee, double _me, double _Mp) {
    Q2 = 4*_Ee*sqrt(Probe.p*Probe.p + _me*_me)*pow(sin(Probe.theta*M_PI/(2*180)), 2);
    nu = _Ee - Probe.p;
    W2 = _Mt*_Mt + 2*_Mt*nu - Q2;
    xB = Q2/(2*_Mp*nu);
}
