/***********************************************
 * Simulation analysis	
 *
 * @author Felix Touchte Codjo
 * @date August 11, 2025
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
    // simu
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/extracted_new_simu.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/test_rec_output2.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/rec_output.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/rec_output2.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/rec_decoded_22092.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/all_simu_rec_output.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    //hipo::bank  trackBank(factory.getSchema("AHDC::track"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nMCtracks =0;
    long unsigned int nKFtracks =0;
    // Histograms
    TH1D* H1_t0 = new TH1D("t0", "t0; time (ns); count", 100, 150, 400);
    TH1D* H1_dt0 = new TH1D("dt0", "leadingEdgeTime - t0; leadingEdgeTime - t0 (ns); count", 100, 0, 400);
    TH1D* H1_time = new TH1D("leadingEdgeTime", "leadingEdgeTime; time (ns); count", 100, 0, 700); 
    TH1D* H1_timeMax = new TH1D("timeMax", "timeMax; timeMax (ns); count", 100, 200, 900); 
    TH1D* H1_deltaTime = new TH1D("deltaTime", "timeMax - leadingEdgeTime (ns); timeMax - leadingEdgeTime (ns); count", 100, 0, 400); 
    TH1D* H1_tot = new TH1D("timeOverThreshold", "timeOverThreshol (ns); timeOverThreshol (ns); count", 100, 150, 752); 
    TH1I* H1_wfType = new TH1I("wfType", "wfType; wfType; count", 6, 0, 6); 
    TH1I* H1_adc = new TH1I("amplitude", "amplitude (adc); count;", 100, 0, 2000); 
    TH2D* H2_times = new TH2D("timeMax, leadingEdgeTime", "timeMax vs leadingEdgeTime; timeMax (ns); leadingEdgeTime (ns);", 10, 200, 900, 100, 0, 700); 
    TH2D* H2_tot_amp = new TH2D("amp, tot", "amplitude vs timeOverThreshold;tot (ns); amp (adc)", 100, 350, 650, 100, 0, 3000); 
    TH2D* H2_deltaTime_adc = new TH2D("deltaTime_adc", "timeMax - leadingEdgeTime vs amplitude;timeMax - leadingEdgeTime (ns); amplitude (ns)", 100, 0, 400, 100, 0, 3700); 
    TH1D* H1_mctime = new TH1D("mctime", "mctime; mctime (ns); count", 100, 0, 400);
    TH1D* H1_diff_mctime = new TH1D("diff_mctime", "mctime - rec_time; mctime - rec_time (ns); count", 100, -200, 200);
    TH1I* H1_nsteps = new TH1I("nsteps", "nsteps in G4; count;", 20, 0, 20);
    std::vector<TH1D*> VecH1_mctime;
    VecH1_mctime.push_back(new TH1D("mctime_nsteps=1", "mctime_nsteps=1; mctime (ns); count", 100, 0, 400));
    VecH1_mctime.push_back(new TH1D("mctime_nsteps=2", "mctime_nsteps=2; mctime (ns); count", 100, 0, 400));
    VecH1_mctime.push_back(new TH1D("mctime_nsteps=3", "mctime_nsteps=3; mctime (ns); count", 100, 0, 400));
    VecH1_mctime.push_back(new TH1D("mctime_nsteps>=4", "mctime_nsteps>=4; mctime (ns); count", 100, 0, 400));
    std::vector<TH1D*> VecH1_tot;
    VecH1_tot.push_back(new TH1D("tot_nsteps=1", "tot_nsteps=1; tot (ns); count", 100, 0, 750));
    VecH1_tot.push_back(new TH1D("tot_nsteps=2", "tot_nsteps=2; tot (ns); count", 100, 0, 750));
    VecH1_tot.push_back(new TH1D("tot_nsteps=3", "tot_nsteps=3; tot (ns); count", 100, 0, 750));
    VecH1_tot.push_back(new TH1D("tot_nsteps>=4", "tot_nsteps>=4; tot (ns); count", 100, 300, 750));
    std::vector<TH1D*> VecH1_amp;
    VecH1_amp.push_back(new TH1D("amp_nsteps=1", "amp_nsteps=1; amp (ns); count", 100, 0, 250));
    VecH1_amp.push_back(new TH1D("amp_nsteps=2", "amp_nsteps=2; amp (ns); count", 100, 0, 350));
    VecH1_amp.push_back(new TH1D("amp_nsteps=3", "amp_nsteps=3; amp (ns); count", 100, 0, 500));
    VecH1_amp.push_back(new TH1D("amp_nsteps>=4", "amp_nsteps>=4; amp (ns); count", 100, 0, 2000));
    std::vector<TH1D*> VecH1_dt0;
    VecH1_dt0.push_back(new TH1D("dt0_nsteps=1", "dt0_nsteps=1; dt0 (ns); count", 100, 0, 400));
    VecH1_dt0.push_back(new TH1D("dt0_nsteps=2", "dt0_nsteps=2; dt0 (ns); count", 100, 0, 400));
    VecH1_dt0.push_back(new TH1D("dt0_nsteps=3", "dt0_nsteps=3; dt0 (ns); count", 100, 0, 400));
    VecH1_dt0.push_back(new TH1D("dt0_nsteps>=4", "dt0_nsteps>=4; dt0 (ns); count", 100, 0, 400));
    TH2D* H2_diffMCtime_amp = new TH2D("corr_diffMCtime_amp", "mctime - rec_time vs amplitude;amp (adc); mctime - rec_time (ns)", 100, 0, 3000, 100, -100, 100); 
    TH2D* H2_diffMCtime_tot = new TH2D("corr_diffMCtime_tot", "mctime - rec_time vs tot;tot (ns); mctime - rec_time (ns)", 100, 400, 700, 100, -100, 100); 
    TH2D* H2_diffMCtime_dt0 = new TH2D("corr_diffMCtime_dt0", "mctime - rec_time vs leadingEdgeTime - t0;leadingEdgeTime - t0 (ns); mctime - rec_time (ns)", 100, 0, 400, 100, -100, 100); 
    // Physics
    TH1D* H1_track_p = new TH1D("track_p", "p; p (MeV); #count", 100, 0, 2000); 
    TH1D* H1_track_pT = new TH1D("track_pT", "pT; pT (MeV); #count", 100, 0, 2000); 
    //TH1D* H1_track_p = new TH1D("track_p", "p; p (MeV); #count", 100, 100, 2000); 
    //TH1D* H1_track_pT = new TH1D("track_pT", "pT; pT (MeV); #count", 100, 100, 2000); 
    TH1D* H1_track_theta = new TH1D("track_theta", "theta; theta (deg); #count", 100, 0, 181); 
    TH1D* H1_track_phi = new TH1D("track_phi", "phi; phi (deg); #count", 100, 0, 361); 
    TH1D* H1_track_vz = new TH1D("track_vz", "vz; vz (cm); #count", 100, -30, 30); 
    //TH1D* H1_track_vz = new TH1D("track_vz", "vz; vz (cm); #count", 100, -50, 15); 
    TH1D* H1_track_sum_adc = new TH1D("track_sum_adc", "#Sigma adc; #Sigma adc; #count", 100, 0, 5000); 
    TH1D* H1_track_res = new TH1D("track_res", "residuals/nhits; residuals/nhits (mm); #count", 100, -1.5, 1); 
    TH2D* H2_track_res_distance = new TH2D("track_res_vs_dist", "residuals vs distance; distance (mm); residuals/nhits (mm)", 100, 0, 4, 100, -1.5, 1); 
    TH1D* H1_track_chi2 = new TH1D("track_chi2", "chi2ndf; chi2ndf (mm^{2}); #count", 100, 0, 8); 
    TH1I* H1_track_nhits = new TH1I("track_nhits", "nhits; nhits; #count", 9, 5, 14); 
    TH1D* H1_hit_distance = new TH1D("hit_distance", "hit distance; distance (mm); #count", 100, 0, 4); 
    TH1D* H1_mc_distance = new TH1D("mc_distance", "distance using mc t2d; distance (mm); #count", 100, 0, 4); // obtained from time 
    TH1I* H1_wfType_bad_distance = new TH1I("wfType_bad_type", "wfType; count;", 6, 0, 6); 
    TH1D* H1_mc_p = new TH1D("mc_p", "p; p (MeV); #count", 100, 230, 270); 
    TH1D* H1_mc_pT = new TH1D("mc_pT", "pT; pT (MeV); #count", 100, 230, 270); 
    TH1D* H1_mc_theta = new TH1D("mc_theta", "theta; theta (deg); #count", 100, 86.38, 87.24); 
    TH1D* H1_mc_phi = new TH1D("mc_phi", "phi; phi (deg); #count", 100, 0, 361); 
    TH1D* H1_mc_vz = new TH1D("mc_vz", "vz; vz (cm); #count", 100, -16, 16); 
    TH2D* H2_pT_amp = new TH2D("corr_pT_amp", "#Sigma_{track} adc vs pT_{mc}; pT_{mc} (MeV); #Sigma_{track} adc", 100, 230, 270, 100, 0, 5000); 
    TH2D* H2_corr_distance = new TH2D("corr_distance", "D_{track} vs D_{mc}; D_{mc} (mm); D_{track} (mm)", 100, 0, 4, 100, 0, 4);
    TH2D* H2_corr_p = new TH2D("corr_p", "p_{track} vs p_{mc}; p_{mc} (MeV); p_{track} (MeV)", 100, 230, 270, 100, 0, 2000);
    TH2D* H2_corr_pT = new TH2D("corr_pT", "pT_{track} vs pT_{mc}; pT_{mc} (MeV); pT_{track} (MeV)", 100, 230, 270, 100, 0, 2000);
    TH2D* H2_corr_phi = new TH2D("corr_phi", "#Phi_{track} vs #Phi_{mc}; #Phi_{mc} (deg); #Phi_{track} (deg)", 100, 0, 361, 100, 0, 361);
    TH2D* H2_corr_theta = new TH2D("corr_theta", "#theta_{track} vs #theta_{mc}; #theta_{mc} (deg); #theta_{track} (deg)", 100, 86.38, 87.24, 100, 0, 181);
    TH2D* H2_corr_vz = new TH2D("corr_vz", "vz_{track} vs vz_{mc}; vz_{mc} (cm); vz_{track} (cm)", 100, -16, 16, 100, -30, 30);
    TH1D* H1_delta_vz = new TH1D("delta_vz", "#Delta vz = vz_{simu} - vz_{track}; #Delta vz (cm); #count", 100, -15, 15); 
    TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta phi = phi_{simu} - phi_{track}; #Delta phi (deg); #count", 100, -30, -5); 


    AhdcCCDB ahdcConstants;
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
        event.getStructure(trackBank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);

        nMCtracks += mcBank.getRows();
        nKFtracks += trackBank.getRows();

        // Loop over waveforms
        for (int i = 0; i < adcBank.getRows(); i++) {
            double sector    = adcBank.getInt("sector", i);
            double layer     = adcBank.getInt("layer", i);
            double component = adcBank.getInt("component", i);
            ahdcT0 obj  = ahdcConstants.get_t0(sector, layer, component);
            double t0 = obj.t0;
            double time = adcBank.getFloat("leadingEdgeTime", i);
            double timeMax = adcBank.getFloat("time", i);
            double tot = adcBank.getFloat("timeOverThreshold", i);
            int wfType = adcBank.getInt("wfType", i);
            int adc = adcBank.getInt("ADC", i);
            //if (wfType == 5) { printf("%d -> %lf ", wfType, time - t0);}
            H2_tot_amp->Fill(tot, adc);
            H1_time->Fill(time);
            H1_t0->Fill(t0);
            H1_dt0->Fill(time - t0);
            H1_timeMax->Fill(timeMax);
            H2_times->Fill(timeMax, time);
            H1_deltaTime->Fill(timeMax-time);
            H1_tot->Fill(tot);
            H1_wfType->Fill(wfType);
            H1_adc->Fill(adc);
            H2_deltaTime_adc->Fill(timeMax-time, adc);
            double mctime = 0.01*wfBank.getInt(std::string("s30").c_str(), i);
            int nsteps = wfBank.getInt(std::string("s29").c_str(), i);
            H1_mctime->Fill(mctime); 
            H1_diff_mctime->Fill(mctime - (time - t0));
            H1_nsteps->Fill(nsteps);
            if (nsteps == 1) {
                VecH1_mctime[0]->Fill(mctime);
                VecH1_tot[0]->Fill(tot);
                VecH1_amp[0]->Fill(adc);
                VecH1_dt0[0]->Fill(time - t0);
            }
            else if (nsteps == 2) {
                VecH1_mctime[1]->Fill(mctime);
                VecH1_tot[1]->Fill(tot);
                VecH1_amp[1]->Fill(adc);
                VecH1_dt0[1]->Fill(time - t0);
            }
            else if (nsteps == 3) {
                VecH1_mctime[2]->Fill(mctime);
                VecH1_tot[2]->Fill(tot);
                VecH1_amp[2]->Fill(adc);
                VecH1_dt0[2]->Fill(time - t0);
            }
            else { // nsteps >= 4
                VecH1_mctime[3]->Fill(mctime);
                VecH1_tot[3]->Fill(tot);
                VecH1_amp[3]->Fill(adc);
                VecH1_dt0[3]->Fill(time - t0);
            }
            H2_diffMCtime_amp->Fill(adc, mctime - (time - t0));
            H2_diffMCtime_tot->Fill(tot, mctime - (time - t0));
            H2_diffMCtime_dt0->Fill(time - t0, mctime - (time - t0));
        }
        // AHDC::kftrack
        // p, pT are already in MeV
        // vz in mm to be converted in cm
        if (trackBank.getRows() > 0) {
            double track_px = trackBank.getFloat("px", 0);
            double track_py = trackBank.getFloat("py", 0);
            double track_pz = trackBank.getFloat("pz", 0);
            double track_vz = trackBank.getFloat("z", 0)*0.1; 
            double track_p, track_theta, track_phi;
            double track_pT;
            futils::cart2polar(track_px, track_py, track_pz, track_p, track_theta, track_phi);
            track_pT = sqrt(pow(track_px,2.0) + pow(track_py, 2.0));
            track_theta = track_theta*180/M_PI;
            track_phi = track_phi*180/M_PI;
            H1_track_p->Fill(track_p); 
            H1_track_pT->Fill(track_pT);
            H1_track_theta->Fill(track_theta);
            H1_track_phi->Fill(track_phi);
            H1_track_vz->Fill(track_vz);
            int nhits = trackBank.getInt("n_hits", 0);
            double res = trackBank.getFloat("sum_residuals", 0);
            double chi2 = trackBank.getFloat("chi2", 0);
            H1_track_nhits->Fill(nhits);
            H1_track_res->Fill(res/nhits);
            H1_track_chi2->Fill(chi2/nhits);
            // for the simu, we know that we have only one track, a rapid hipo-utils...
            for (int hit_row = 0; hit_row < hitBank.getRows(); hit_row++) {
                double hit_distance = hitBank.getDouble("doca", hit_row);
                double hit_time = hitBank.getDouble("time", hit_row);
                int wfType = adcBank.getInt("wfType", hitBank.getInt("id", hit_row) - 1);
                double mctime = 0.01*wfBank.getInt(std::string("s30").c_str(), hitBank.getInt("id", hit_row) - 1);
                if (hit_distance < 0) {
                    //printf("%d --> %lf ", wfType, hit_distance);
                    H1_wfType_bad_distance->Fill(wfType);
                }
                H1_hit_distance->Fill(hit_time < 0 ? 0 : hit_distance);
                double mc_distance = (mctime < 0) ? 0 : -0.0497 -0.00667*mctime + 0.389*sqrt(mctime) - 0.189*pow(mctime, 1.0/3);
                H1_mc_distance->Fill(mc_distance);
                H2_corr_distance->Fill(mc_distance, hit_distance);
                H2_track_res_distance->Fill(hit_distance, res/nhits);
            }
            // MC::Particle or AHDC::mc
            // p, pT to be converted in MeV
            // vz already in cm
            double mc_px = 1000*mcBank.getFloat("px", 0);
            double mc_py = 1000*mcBank.getFloat("py", 0);
            double mc_pz = 1000*mcBank.getFloat("pz", 0);
            double mc_vz = mcBank.getFloat("vz", 0); 
            double mc_p, mc_theta, mc_phi;
            double mc_pT;
            futils::cart2polar(mc_px, mc_py, mc_pz, mc_p, mc_theta, mc_phi);
            mc_pT = sqrt(pow(mc_px,2.0) + pow(mc_py, 2.0));
            mc_theta = mc_theta*180/M_PI;
            mc_phi = mc_phi*180/M_PI;
            H1_mc_p->Fill(mc_p); 
            H1_mc_pT->Fill(mc_pT);
            H1_mc_theta->Fill(mc_theta);
            H1_mc_phi->Fill(mc_phi);
            H1_mc_vz->Fill(mc_vz);
            // Correlation
            H1_delta_vz->Fill(mc_vz - track_vz);
            H1_delta_phi->Fill(mc_phi - track_phi);
            H2_pT_amp->Fill(mc_pT, trackBank.getInt("sum_adc", 0));
            H1_track_sum_adc->Fill(trackBank.getInt("sum_adc", 0));
            H2_corr_phi->Fill(mc_phi, track_phi);
            H2_corr_theta->Fill(mc_theta, track_theta);
            H2_corr_p->Fill(mc_p, track_p);
            H2_corr_pT->Fill(mc_pT, track_pT);
            H2_corr_vz->Fill(mc_vz, track_vz);
        }
    }
    
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    const char * output = "./output/new_simu_data.root";
    TFile *f = new TFile(output, "RECREATE");
    
    TDirectory *time_study = f->mkdir("time_study");
    time_study->cd();
// >>>>>>>>>   Process H2_deltaTime_adc
    TCanvas* canvas1 = new TCanvas("c1_deltaTime_amplitude", "deltaTime vs amplitude; timeMax - leadingedgeTime (ns); amplitude (adc)");
    canvas1->SetLogz();
    H2_deltaTime_adc->Draw("colz");
    TText* text1 = new TText();
    text1->SetTextSize(0.02);
    text1->SetTextColor(kRed);
    int nbins = H2_deltaTime_adc->GetXaxis()->GetNbins();
    printf("> Analyse correlation between deltaTime (ns) and amplitude (adc)\n");
    printf(">   nbins x axis : %d\n", nbins);
    int proj_nbins = 0.05*nbins;
    int Npts = nbins/proj_nbins;
    TGraphErrors* gr = new TGraphErrors(Npts);
    for (int i = 0; i < Npts; i++) {
        TH1D* H1_projTime = H2_deltaTime_adc->ProjectionX(TString::Format("_px_%d", i).Data(), i*proj_nbins, (i+1)*proj_nbins);
        //H1_projTime->Write(TString::Format("_px_%d", i).Data());
        double mean = H1_projTime->GetMean();
        double stdDev = H1_projTime->GetStdDev();
        printf(">   mean : %lf, stdDev : %lf\n", mean, stdDev);
        double y_bin_inf = H2_deltaTime_adc->GetYaxis()->GetBinCenter(i*proj_nbins); // adc_inf for this collection of bins
        double y_bin_sup = H2_deltaTime_adc->GetYaxis()->GetBinCenter((i+1)*proj_nbins); // adc_sup this collection of bins
        gr->SetPoint(i, mean, 0.5*(y_bin_inf+y_bin_sup));
        gr->SetPointError(i, stdDev, 0);
        text1->DrawText(mean + stdDev, 0.5*(y_bin_inf+y_bin_sup), TString::Format("%.2lf +/- %.2lf", mean, stdDev).Data());
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
    H1_time->Write("leadingEdgeTime");
    H1_t0->Write("t0");
    H1_timeMax->Write("timeMax");
    H1_deltaTime->Write("dt");
    H2_times->Write("corr_timeMax_leadingEdgeTime");
    H2_tot_amp->Write("corr_tot_amp");
    H1_wfType->Write("wfType");
    H2_deltaTime_adc->Write("corr_dt_amp");
    H1_mctime->Write("mctime");
    H1_diff_mctime->Write("diff_mctime");
    H1_nsteps->Write("nsteps_G4");
    VecH1_mctime[0]->Write("mctime_nsteps=1");
    VecH1_mctime[1]->Write("mctime_nsteps=2");
    VecH1_mctime[2]->Write("mctime_nsteps=3");
    VecH1_mctime[3]->Write("mctime_nsteps>=4");
    H1_tot->Write("tot");
    VecH1_tot[0]->Write("tot_nsteps=1");
    VecH1_tot[1]->Write("tot_nsteps=2");
    VecH1_tot[2]->Write("tot_nsteps=3");
    VecH1_tot[3]->Write("tot_nsteps>=4");
    H1_adc->Write("amp");
    VecH1_amp[0]->Write("amp_nsteps=1");
    VecH1_amp[1]->Write("amp_nsteps=2");
    VecH1_amp[2]->Write("amp_nsteps=3");
    VecH1_amp[3]->Write("amp_nsteps>=4");
    H1_dt0->Write("dt0");
    VecH1_dt0[0]->Write("dt0_nsteps=1");
    VecH1_dt0[1]->Write("dt0_nsteps=2");
    VecH1_dt0[2]->Write("dt0_nsteps=3");
    VecH1_dt0[3]->Write("dt0_nsteps>=4");
    H2_diffMCtime_amp->Write("corr_diffMCtime_amp");
    H2_diffMCtime_tot->Write("corr_diffMCtime_tot");
    H2_diffMCtime_dt0->Write("corr_diffMCtime_dt0");

    // others
    TDirectory *track_study = f->mkdir("track_study");
    track_study->cd();
    H1_track_p->Write("track_p");
    H1_track_pT->Write("track_pT");
    H1_track_theta->Write("track_theta");
    H1_track_phi->Write("track_phi");
    H1_track_vz->Write("track_vz");
    H1_track_nhits->Write("nhits");
    H1_track_res->Write("residuals");
    H1_track_chi2->Write("chi2ndf");
    H1_hit_distance->Write("recon_hit_distance");
    H1_mc_distance->Write("recon_mc_distance");
    H2_corr_distance->Write("corr_distance");
    H2_track_res_distance->Write("corr_residuals_distance");
    H1_wfType_bad_distance->Write("wfType_bad_distance");
    H1_mc_p->Write("mc_p");
    H1_mc_pT->Write("mc_pT");
    H1_mc_theta->Write("mc_theta");
    H1_mc_phi->Write("mc_phi");
    H1_mc_vz->Write("mc_vz");
// >>>>>>>>>   Compare track and simu
    TCanvas* canvas2 = new TCanvas("c1_track_study", "Simulation vs Reconstruction");
    canvas2->Divide(3,2);
    canvas2->cd(1);
    THStack* stack1 = new THStack("stack_p", "#bf{p} #color[4]{Reconstruction} vs #color[2]{Simulation}; p (MeV)");
    stack1->Add(H1_track_p);
    H1_mc_p->SetLineColor(kRed);
    stack1->Add(H1_mc_p);
    stack1->Draw("nostack");
    canvas2->cd(2);
    THStack* stack2 = new THStack("stack_pT", "#bf{pT} #color[4]{Reconstruction} vs #color[2]{Simulation}; pT (MeV)");
    stack2->Add(H1_track_pT);
    H1_mc_pT->SetLineColor(kRed);
    stack2->Add(H1_mc_pT);
    stack2->Draw("nostack");
    canvas2->cd(3);
    THStack* stack3 = new THStack("stack_theta", "#bf{theta} #color[4]{Reconstruction} vs #color[2]{Simulation}; theta (deg)");
    stack3->Add(H1_track_theta);
    H1_mc_theta->SetLineColor(kRed);
    stack3->Add(H1_mc_theta);
    stack3->Draw("nostack");
    canvas2->cd(4);
    THStack* stack4 = new THStack("stack_phi", "#bf{phi} #color[4]{Reconstruction} vs #color[2]{Simulation}; phi (deg)");
    stack4->Add(H1_track_phi);
    H1_mc_phi->SetLineColor(kRed);
    stack4->Add(H1_mc_phi);
    stack4->Draw("nostack");
    canvas2->cd(5);
    THStack* stack5 = new THStack("stack_vz", "#bf{vz} #color[4]{Reconstruction} vs #color[2]{Simulation}; vz (cm)");
    stack5->Add(H1_track_vz);
    H1_mc_vz->SetLineColor(kRed);
    stack5->Add(H1_mc_vz);
    stack5->Draw("nostack");
    canvas2->Write("Side_by_side_comparaison");
// <<<<<<<<<<<<<< Compare track and simu
    H2_pT_amp->Write("corr_pT_ADC");
    H1_track_sum_adc->Write("track_sum_adc");
    H2_corr_p->Write("corr_p");
    H2_corr_pT->Write("corr_pT");
    H2_corr_theta->Write("corr_theta");
    H2_corr_phi->Write("corr_phi");
    H2_corr_vz->Write("corr_vz");
    H1_delta_phi->Write("delta_phi");
    H1_delta_vz->Write("delta_vz");

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
