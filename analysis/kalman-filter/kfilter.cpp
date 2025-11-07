/***********************************************
 * Study the Kalman Filter for the AHDC
 *
 * @author Felix Touchte Codjo
 * @date October 23, 2025
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
    
    //const char * output = "./output/kfilter-v1.root";
    //const char * output = "./output/kfilter-deuteron-v2.root";
    //const char * output = "./output/kfilter-deuteron-v4.root";
    const char * output = "./output/kfilter-deuteron-v17.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-proton-v1.hipo";
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v17.hipo";
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
    // adc
    TH1D* H1_occ = new TH1D("occupancy", "occupancy; wire; count [%]", 576, 0, 576);
    TH1D* H1_t0 = new TH1D("t0", "t0; time (ns); count", 100, 150, 400);
    TH1D* H1_leadingEdgeTime = new TH1D("leadingEdgeTime", "leadingEdgeTime; leadingEdgeTime (ns); count", 100, 100, 600); 
    TH1D* H1_tot = new TH1D("timeOverThreshold", "timeOverThreshol (ns); timeOverThreshol (ns); count", 100, 150, 750); 
    TH1I* H1_wfType = new TH1I("wfType", "wfType; wfType; count", 6, 0, 6); 
    TH1I* H1_adc = new TH1I("amplitude", "amplitude (adc); count;", 100, 0, 3000); 
    TH2D* H2_tot_amp = new TH2D("amp, tot", "amplitude vs timeOverThreshold;tot (ns); amp (adc)", 100, 150, 750, 100, 0, 3000); 
    TH1I* H1_nsteps = new TH1I("nsteps", "nsteps in G4; nsteps; count;", 20, 0, 20);
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
    // track vs mc
    TH1D* H1_track_p = new TH1D("track_p", "p; p (MeV); #count", 100, 0, 1000); 
    TH1D* H1_track_pT = new TH1D("track_pT", "pT; pT (MeV); #count", 100, 0, 1000); 
    TH1D* H1_track_theta = new TH1D("track_theta", "theta; theta (deg); #count", 100, 0, 181); 
    TH1D* H1_track_phi = new TH1D("track_phi", "phi; phi (deg); #count", 100, 0, 361); 
    TH1D* H1_track_vz = new TH1D("track_vz", "vz; vz (cm); #count", 100, -16, 16); 
    TH1D* H1_track_sum_adc = new TH1D("track_sum_adc", "#Sigma adc; #Sigma adc; #count", 100, 0, 5000); 
    TH1D* H1_track_sum_residuals = new TH1D("track_res", "Sum residuals; Sum residuals (mm); #count", 100, -5, 5); 
    TH1D* H1_track_chi2 = new TH1D("track_chi2", "chi2ndf; chi2ndf (mm^{2}); #count", 100, 0, 3); 
    TH1I* H1_track_nhits = new TH1I("track_nhits", "nhits; nhits; #count", 9, 5, 14); 
    TH2D* H2_track_pT_amp = new TH2D("corr_pT_amp", "#Sigma_{track} adc vs pT_{track}; pT_{track} (MeV); #Sigma_{track} adc", 100, 225, 275, 100, 0, 5000); 
    TH1D* H1_mc_p = new TH1D("mc_p", "p; p (MeV); #count", 100, 0, 1000); 
    TH1D* H1_mc_pT = new TH1D("mc_pT", "pT; pT (MeV); #count", 100, 0, 1000); 
    TH1D* H1_mc_theta = new TH1D("mc_theta", "theta; theta (deg); #count", 100, 0, 181); 
    TH1D* H1_mc_phi = new TH1D("mc_phi", "phi; phi (deg); #count", 100, 0, 361); 
    TH1D* H1_mc_vz = new TH1D("mc_vz", "vz; vz (cm); #count", 100, -16, 16); 
    TH2D* H2_corr_p = new TH2D("corr_p", "p_{track} vs p_{mc}; p_{mc} (MeV); p_{track} (MeV)", 100, 190, 310, 100, 190, 310);
    TH2D* H2_corr_pT = new TH2D("corr_pT", "pT_{track} vs pT_{mc}; pT_{mc} (MeV); pT_{track} (MeV)", 100, 190, 310, 100, 190, 310);
    TH2D* H2_corr_phi = new TH2D("corr_phi", "#Phi_{track} vs #Phi_{mc}; #Phi_{mc} (deg); #Phi_{track} (deg)", 100, 0, 361, 100, 0, 361);
    TH2D* H2_corr_theta = new TH2D("corr_theta", "#theta_{track} vs #theta_{mc}; #theta_{mc} (deg); #theta_{track} (deg)", 100, 55, 125, 100, 55, 125);
    TH2D* H2_corr_vz = new TH2D("corr_vz", "vz_{track} vs vz_{mc}; vz_{mc} (cm); vz_{track} (cm)", 100, -16, 16, 100, -15, 15);
    TH1D* H1_delta_vz = new TH1D("delta_vz", "#Delta vz = vz_{mc} - vz_{track}; #Delta vz (cm); #count", 100, -5, 10); 
    TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta phi = phi_{mc} - phi_{track}; #Delta phi (deg); #count", 100, -5, 5); 
    // hit 
    TH1D* H1_distance = new TH1D("distance", "hit distance; distance (mm); #count", 100, 0, 4); 
    TH1D* H1_mcdistance = new TH1D("mcdistance", "distance using mc t2d; distance (mm); #count", 100, 0, 4); // obtained from time 
    TH1D* H1_deltaDistance = new TH1D("deltaDistance", "deltaDistance = mcdistance - distance; deltaDistance (mm); #count", 100, -1.5, 1.5); 
    TH1D* H1_time = new TH1D("time", "time = leadingEdgeTime - t0 - startTime (= 0 for simu); time (ns); count", 100, 0, 250);
    TH1D* H1_mctime = new TH1D("mctime", "mctime; mctime (ns); count", 100, 0, 250);
    TH1D* H1_deltaTime = new TH1D("deltaTime", "mctime - rectime; mctime - rectime (ns); count", 100, -40, 60);
    TH2D* H2_deltaTime_amp = new TH2D("corr_deltaTime_amp", "deltaTime vs amplitude; mctime - rectime = deltaTime (ns); amp (adc)", 100, -40, 60, 100, 0, 3000); 
    TH2D* H2_deltaTime_tot = new TH2D("corr_deltaTime_tot", "deltaTime vs timeOverThreshold; mctime - rectime = deltaTime (ns); tot (ns)", 100, -40, 60, 100, 150, 750); 
    TH2D* H2_corr_distance = new TH2D("corr_distance", "distance vs mcdistance; mcdistance (mm); distance (mm)", 100, 0, 4, 100, 0, 4);
    TH2D* H2_corr_time = new TH2D("corr_time", "time vs mctime; mctime (mm); time (mm)", 100, 0, 250, 100, 0, 250);
    TH1D* H1_residual = new TH1D("residuals", "residuals; residuals (mm); #count", 100, -4, 4); 
    TH2D* H2_time2distance = new TH2D("corr_time2distance", "time2distance; time (ns); KF distance (mm)", 100, 0, 250, 100, 0, 4); 
    std::vector<TH1D*> VecH1_noise;
    VecH1_noise.push_back(new TH1D("noise s0", "s0; s0 (adc); count", 100, 250, 380));
    VecH1_noise.push_back(new TH1D("noise s1", "s1; s1 (adc); count", 100, 250, 380));
    VecH1_noise.push_back(new TH1D("noise s2", "s2; s2 (adc); count", 100, 250, 380));
    VecH1_noise.push_back(new TH1D("noise s3", "s3; s3 (adc); count", 100, 250, 380));

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
            H1_track_nhits->Fill(nhits);
            H1_track_sum_residuals->Fill(trackBank.getFloat("sum_residuals", 0));
            H1_track_chi2->Fill(trackBank.getFloat("chi2", 0)/nhits);
            int trackid = trackBank.getInt("trackid", 0);
            //printf("trackid : %d \n", trackid);
            // for the simu, we know that we have only one track, a rapid hipo-utils...
            for (int row = 0; row < hitBank.getRows(); row++) {
                if (trackid != hitBank.getInt("trackid", row)) continue;
                int adcId = hitBank.getInt("id", row) - 1;
                int    wfType = adcBank.getInt("wfType", adcId);
                double mctime = 0.01*wfBank.getInt(std::string("s30").c_str(), adcId);
                double mcdistance = (mctime < 0) ? 0 : -0.0497 -0.00667*mctime + 0.389*sqrt(mctime) - 0.189*pow(mctime, 1.0/3);
                double time = hitBank.getDouble("time", row);
                double distance = hitBank.getDouble("doca", row);
                H1_distance->Fill(distance);
                H1_mcdistance->Fill(mcdistance);
                H1_time->Fill(time);
                H1_mctime->Fill(mctime);
                H1_deltaDistance->Fill(mcdistance-distance);
                H1_deltaTime->Fill(mctime-time);
                H2_corr_distance->Fill(mcdistance, distance);
                H2_corr_time->Fill(mctime, time);
                double residual = hitBank.getDouble("residual", row);
                H1_residual->Fill(residual);
                // residual = hit.doca - kf.doca
                H2_time2distance->Fill(time, distance - residual);
                ///////////////////
                // AHDC::adc
                // ////////////////
                //printf("hit row : %d , adc id : %d\n", row, adcId);
                double sector    = adcBank.getInt("sector", adcId);
                double layer     = adcBank.getInt("layer", adcId);
                double component = adcBank.getInt("component", adcId);
                ahdcT0 obj  = ahdcConstants.get_t0(sector, layer, component);
                double t0 = obj.t0;
                double leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", adcId);
                double tot             = adcBank.getFloat("timeOverThreshold", adcId);
                int    adc             = adcBank.getInt("ADC", adcId);
                int    nsteps = wfBank.getInt(std::string("s29").c_str(), adcId);
                // Histograms
                H1_occ->Fill(AhdcCCDB::wireUniqueId(sector, layer, component));
                H1_leadingEdgeTime->Fill(leadingEdgeTime);
                H1_t0->Fill(t0);
                H1_tot->Fill(tot);
                H1_wfType->Fill(wfType);
                H1_adc->Fill(adc);
                H1_nsteps->Fill(nsteps);
                H2_tot_amp->Fill(tot, adc);
                H2_deltaTime_amp->Fill(mctime - time, adc);
                H2_deltaTime_tot->Fill(mctime - time, tot);
                    // Noise level
                int s0 = wfBank.getShort("s1", adcId);
                int s1 = wfBank.getShort("s2", adcId);
                int s2 = wfBank.getShort("s3", adcId);
                int s3 = wfBank.getShort("s4", adcId);
                VecH1_noise[0]->Fill(s0);
                VecH1_noise[1]->Fill(s1);
                VecH1_noise[2]->Fill(s2);
                VecH1_noise[3]->Fill(s3);
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
            H2_track_pT_amp->Fill(track_pT, trackBank.getInt("sum_adc", 0));
            H1_track_sum_adc->Fill(trackBank.getInt("sum_adc", 0));
            H2_corr_phi->Fill(mc_phi, track_phi);
            H2_corr_theta->Fill(mc_theta, track_theta);
            H2_corr_p->Fill(mc_p, track_p);
            H2_corr_pT->Fill(mc_pT, track_pT);
            H2_corr_vz->Fill(mc_vz, track_vz);
        }
    }
    H1_occ->Scale(100.0/nKFtracks); 
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    TFile *f = new TFile(output, "RECREATE");
    
    TDirectory *adc_dir = f->mkdir("adc");
    adc_dir->cd();
    H1_occ->Write("occupancy");
    H1_wfType->Write("wfType");
    H1_leadingEdgeTime->Write("leadingEdgeTime");
    H1_tot->Write("timeOverThreshold");
    //VecH1_tot[0]->Write("tot_nsteps=1");
    //VecH1_tot[1]->Write("tot_nsteps=2");
    //VecH1_tot[2]->Write("tot_nsteps=3");
    //VecH1_tot[3]->Write("tot_nsteps>=4");
    H1_adc->Write("amplitude");
    //VecH1_amp[0]->Write("amp_nsteps=1");
    //VecH1_amp[1]->Write("amp_nsteps=2");
    //VecH1_amp[2]->Write("amp_nsteps=3");
    //VecH1_amp[3]->Write("amp_nsteps>=4");
    H2_tot_amp->Write("corr_tot_amp");
    H1_nsteps->Write("nsteps_G4");
    H1_t0->Write("t0");
    
    // hit
    TDirectory *hit_dir = f->mkdir("hit");
    hit_dir->cd();
    H1_residual->Write("residual");
    H2_time2distance->Write("time2distance");
    H1_distance->Write("distance");
    H1_mcdistance->Write("mcdistance");
    H1_deltaDistance->Write("delta_distance");
    H1_time->Write("time");
    H1_mctime->Write("mctime");
    H1_deltaTime->Write("delta_time");
    H2_corr_distance->Write("corr_distance");
    H2_corr_time->Write("corr_time");
    H2_deltaTime_amp->Write("corr_deltaTime_amp");
    H2_deltaTime_tot->Write("corr_deltaTime_tot");
    //VecH1_mctime[0]->Write("mctime_nsteps=1");
    //VecH1_mctime[1]->Write("mctime_nsteps=2");
    //VecH1_mctime[2]->Write("mctime_nsteps=3");
    //VecH1_mctime[3]->Write("mctime_nsteps>=4");
    
    // track
    TDirectory *track_dir = f->mkdir("track");
    track_dir->cd();
    H1_track_p->Write("track_p");
    H1_track_pT->Write("track_pT");
    H1_track_theta->Write("track_theta");
    H1_track_phi->Write("track_phi");
    H1_track_vz->Write("track_vz");
    H1_track_nhits->Write("nhits");
    H1_track_sum_residuals->Write("sum_residuals");
    H1_track_chi2->Write("chi2ndf");
    H1_track_sum_adc->Write("track_sum_adc");
    H2_track_pT_amp->Write("corr_pT_ADC");
    // mctrack
    TDirectory *mctrack_dir = f->mkdir("mctrack");
    mctrack_dir->cd();
    H1_mc_p->Write("mc_p");
    H1_mc_pT->Write("mc_pT");
    H1_mc_theta->Write("mc_theta");
    H1_mc_phi->Write("mc_phi");
    H1_mc_vz->Write("mc_vz");
    // track vs mc
    TDirectory *comp_dir = f->mkdir("comparison_track");
    comp_dir->cd();
    {
    // Compare track and simu
        TCanvas* canvas2 = new TCanvas("c1_track_study", "Simulated vs Reconstructed");
        canvas2->Divide(3,2);
        canvas2->cd(1);
        THStack* stack1 = new THStack("stack_p", "#bf{p} #color[4]{Reconstructed} vs #color[2]{Simulated}; p (MeV)");
        stack1->Add(H1_track_p);
        H1_mc_p->SetLineColor(kRed);
        stack1->Add(H1_mc_p);
        stack1->Draw("nostack");
        canvas2->cd(2);
        THStack* stack2 = new THStack("stack_pT", "#bf{pT} #color[4]{Reconstructed} vs #color[2]{Simulated}; pT (MeV)");
        stack2->Add(H1_track_pT);
        H1_mc_pT->SetLineColor(kRed);
        stack2->Add(H1_mc_pT);
        stack2->Draw("nostack");
        canvas2->cd(3);
        THStack* stack3 = new THStack("stack_theta", "#bf{theta} #color[4]{Reconstructed} vs #color[2]{Simulated}; theta (deg)");
        stack3->Add(H1_track_theta);
        H1_mc_theta->SetLineColor(kRed);
        stack3->Add(H1_mc_theta);
        stack3->Draw("nostack");
        canvas2->cd(4);
        THStack* stack4 = new THStack("stack_phi", "#bf{phi} #color[4]{Reconstructed} vs #color[2]{Simulated}; phi (deg)");
        stack4->Add(H1_track_phi);
        H1_mc_phi->SetLineColor(kRed);
        stack4->Add(H1_mc_phi);
        stack4->Draw("nostack");
        canvas2->cd(5);
        THStack* stack5 = new THStack("stack_vz", "#bf{vz} #color[4]{Reconstructed} vs #color[2]{Simulated}; vz (cm)");
        stack5->Add(H1_track_vz);
        H1_mc_vz->SetLineColor(kRed);
        stack5->Add(H1_mc_vz);
        stack5->Draw("nostack");
        canvas2->Write("Side_by_side_comparaison");
    }
    H2_corr_p->Write("corr_p");
    H2_corr_pT->Write("corr_pT");
    H2_corr_theta->Write("corr_theta");
    H2_corr_phi->Write("corr_phi");
    H2_corr_vz->Write("corr_vz");
    H1_delta_phi->Write("delta_phi");
    H1_delta_vz->Write("delta_vz");
    // noise
    TDirectory *noise_dir = f->mkdir("noise_dir");
    noise_dir->cd();
    VecH1_noise[0]->Write("s0");
    VecH1_noise[1]->Write("s1");
    VecH1_noise[2]->Write("s2");
    VecH1_noise[3]->Write("s3");


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
