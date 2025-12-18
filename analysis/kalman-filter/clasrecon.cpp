/***********************************************
 * Study CLAS recon reconstruction
 *
 * @author Felix Touchte Codjo
 * @date November 27, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>

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
#include "TMultiGraph.h"

#include "futils.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities
void computeSphericalVariance(double mu_x, double mu_y, double mu_z, double var_x, double var_y, double var_z, double & var_r, double & var_theta, double & var_phi, const char * nature = "");


int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
   const char * output = "./output/clas-recon-v36.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v36.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  kfmonBank(factory.getSchema("AHDC::kftrack:mon"));
    hipo::bank  track0Bank(factory.getSchema("AHDC::track"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::bank  trueBank(factory.getSchema("MC::True"));
    hipo::bank  runConfigBank(factory.getSchema("RUN::config"));
    hipo::bank  recBank(factory.getSchema("REC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nMCtracks =0;
    long unsigned int nKFtracks =0;
   

    // Histograms
    TH2D* H2_corr_phi = new TH2D("corr_phi", "#phi_{mc} vs #phi_{rec}; #phi_{rec} (deg); #phi_{mc} (deg)", 100, 0, 361, 100, 0, 361); 
    TH2D* H2_corr_theta = new TH2D("corr_theta", "#theta_{mc} vs #theta_{rec}; #theta_{rec} (deg); #theta_{mc} (deg)", 100, 4, 10, 100, 4, 10); 
    TH2D* H2_corr_p = new TH2D("corr_p", "p_{mc} vs p_{rec}; p_{rec} (GeV); p_{mc} (GeV)", 100, 1.8, 2.5, 100, 1.8, 2.5); 
    TH2D* H2_corr_vz = new TH2D("corr_vz", "vz_{mc} vs vz_{rec}; vz_{rec} (cm); vz_{mc} (cm)", 100, -16, 16, 100, -16, 16); 
    TH1D* H1_mc_p = new TH1D("mc_p", "p; p (GeV); #count", 100, 2.15, 2.2);
    TH1D* H1_mc_theta = new TH1D("mc_theta", "theta; theta (deg); #count", 100, 4, 10);
    TH1D* H1_mc_phi = new TH1D("mc_phi", "phi; phi (deg); #count", 100, 0, 361);
    TH1D* H1_mc_vx = new TH1D("mc_vx", "vx; vx (cm); #count", 100, -16, 16);
    TH1D* H1_mc_vy = new TH1D("mc_vy", "vy; vy (cm); #count", 100, -16, 16);
    TH1D* H1_mc_vz = new TH1D("mc_vz", "vz; vz (cm); #count", 100, -16, 16);
    TH1D* H1_rec_p = new TH1D("rec_p", "p; p (GeV); #count", 100, 2.15, 2.2);
    TH1D* H1_rec_theta = new TH1D("rec_theta", "theta; theta (deg); #count", 100, 4, 10);
    TH1D* H1_rec_phi = new TH1D("rec_phi", "phi; phi (deg); #count", 100, 0, 361);
    TH1D* H1_rec_vx = new TH1D("rec_vx", "vx; vx (cm); #count", 100, -16, 16);
    TH1D* H1_rec_vy = new TH1D("rec_vy", "vy; vy (cm); #count", 100, -16, 16);
    TH1D* H1_rec_vz = new TH1D("rec_vz", "vz; vz (cm); #count", 100, -16, 16);
    TH1D* H1_delta_vz = new TH1D("delta_vz", "#Delta vz = vz_{mc} - vz_{rec}; #Delta vz (cm); #count", 100, -16, 16);
    TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta #phi = #phi_{mc} - #phi_{rec}; #Delta #phi (deg); #count", 100, -10, 10);
    TH1D* H1_delta_theta = new TH1D("delta_theta", "#Delta #theta = #theta_{mc} - #theta_{rec}; #Delta #theta (deg); #count", 100, -2, 2);

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
        event.getStructure(kfmonBank);
        event.getStructure(track0Bank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);
        event.getStructure(trueBank);
        event.getStructure(runConfigBank);
        event.getStructure(recBank);
        
        //if (recBank.getRows() < 1) continue; 

        nMCtracks += mcBank.getRows();
        nKFtracks += trackBank.getRows();

        // Loop over REC::Particle
        for (int i = 0; i < recBank.getRows(); i++) {
            if ((int) recBank.get("pid", i) == 11) {
                double rec_px = recBank.get("px", i); 
                double rec_py = recBank.get("py", i); 
                double rec_pz = recBank.get("pz", i); 
                double rec_p; 
                double rec_theta; 
                double rec_phi;
                futils::cart2polarDEG(rec_px, rec_py, rec_pz, rec_p, rec_theta, rec_phi);
                double rec_vx = recBank.get("vx", i); 
                double rec_vy = recBank.get("vy", i); 
                double rec_vz = recBank.get("vz", i); 
                // Fill rec histograms
                H1_rec_p->Fill(rec_p); 
                H1_rec_theta->Fill(rec_theta); 
                H1_rec_phi->Fill(rec_phi); 
                H1_rec_vz->Fill(rec_vz); 
                H1_rec_vy->Fill(rec_vy); 
                H1_rec_vx->Fill(rec_vx); 
                // Loop over MC::Particle
                for (int j = 0; j < mcBank.getRows(); j++) {
                    if ((int) mcBank.get("pid", i) == 11) {
                        double mc_px = mcBank.get("px", i); 
                        double mc_py = mcBank.get("py", i); 
                        double mc_pz = mcBank.get("pz", i); 
                        double mc_p; 
                        double mc_theta; 
                        double mc_phi;
                        futils::cart2polarDEG(mc_px, mc_py, mc_pz, mc_p, mc_theta, mc_phi);
                        double mc_vx = mcBank.get("vx", i); 
                        double mc_vy = mcBank.get("vy", i); 
                        double mc_vz = mcBank.get("vz", i); 
                        // Fill mc histogram
                        H1_mc_p->Fill(mc_p); 
                        H1_mc_theta->Fill(mc_theta); 
                        H1_mc_phi->Fill(mc_phi); 
                        H1_mc_vz->Fill(mc_vz); 
                        H1_mc_vy->Fill(mc_vy); 
                        H1_mc_vx->Fill(mc_vx); 
                        // Fill Coorelations
                        H2_corr_p->Fill(rec_p, mc_p);
                        H2_corr_theta->Fill(rec_theta, mc_theta);
                        H2_corr_phi->Fill(rec_phi, mc_phi);
                        H2_corr_vz->Fill(rec_vz, mc_vz);
                        H1_delta_vz->Fill(mc_vz-rec_vz);
                        H1_delta_theta->Fill(mc_theta-rec_theta);
                        H1_delta_phi->Fill(mc_phi-rec_phi);
                        break; // first deuteron
                    }
                } // end loop over mc particle
                break; // first electron
            }
        } // end loop over rec particle
        
    } // end loop over events
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    TFile *f = new TFile(output, "RECREATE");
    TDirectory * mc_dir = f->mkdir("mc");  
    TDirectory * rec_dir = f->mkdir("rec");  
    TDirectory * corr_dir = f->mkdir("corr");  
    // mc
    mc_dir->cd();
    H1_mc_p->Write(H1_mc_p->GetName());
    H1_mc_theta->Write(H1_mc_theta->GetName());
    H1_mc_phi->Write(H1_mc_phi->GetName());
    H1_mc_vz->Write(H1_mc_vz->GetName());
    H1_mc_vy->Write(H1_mc_vy->GetName());
    H1_mc_vx->Write(H1_mc_vx->GetName());
    // rec
    rec_dir->cd();
    H1_rec_p->Write(H1_rec_p->GetName());
    H1_rec_theta->Write(H1_rec_theta->GetName());
    H1_rec_phi->Write(H1_rec_phi->GetName());
    H1_rec_vz->Write(H1_rec_vz->GetName());
    H1_rec_vy->Write(H1_rec_vy->GetName());
    H1_rec_vx->Write(H1_rec_vx->GetName());
    // corr
    corr_dir->cd();
    H2_corr_phi->Write("corr_phi");
    H2_corr_theta->Write("corr_theta");
    H2_corr_p->Write("corr_p");
    H2_corr_vz->Write("corr_vz");
    H1_delta_vz->Write(H1_delta_vz->GetName());
    H1_delta_theta->Write(H1_delta_theta->GetName());
    H1_delta_phi->Write(H1_delta_phi->GetName());

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

void computeSphericalVariance(double mu_x, double mu_y, double mu_z, double var_x, double var_y, double var_z, double & var_r, double & var_theta, double & var_phi, const char * nature) {
    // prevent NaN ot Inf values
    // if mu_x and mu_y are strictly equal to 0 (at beginning of the Kalman Filter for example)
    // phi can not be defined
    // theta is equals to 0 but its derivative is infinite at 0
    // we then assume that all the variance if given by the z component to the r component
    // I assume x, y and z independents but in reality we can correlations
    if ((std::abs(mu_x) < 1e-9) && (std::abs(mu_y) < 1e-9)) {
        var_r = var_z;
        var_theta = 0;
        var_phi = 0;
    } else {
        // var_r
        double r = sqrt(mu_x*mu_x+mu_y*mu_y+mu_z*mu_z);
        double drdx = mu_x/r;
        double drdy = mu_y/r;
        double drdz = mu_z/r;
        var_r = drdx*drdx*var_x + drdy*drdy*var_y + drdz*drdz*var_z;
        // var_theta
        double dthetadz = (-1/r)/sqrt(1-pow(mu_z/r,2.0));
        double dthetady = (1/r)*(mu_z/r)*(mu_y/r)/sqrt(1-pow(mu_z/r,2.0));
        double dthetadx = (1/r)*(mu_z/r)*(mu_x/r)/sqrt(1-pow(mu_z/r,2.0));
        var_theta = dthetadx*dthetadx*var_x + dthetady*dthetady*var_y + dthetadz*dthetadz*var_z;
        var_theta *= pow(180.0/M_PI,2.0); // convert to deg²
        // var_phi
        double dphidz = 0;
        double rho2 = mu_x*mu_x+mu_y*mu_y;
        double dphidy = mu_x/rho2;
        double dphidx = -mu_y/rho2;
        var_phi = dphidx*dphidx*var_x + dphidy*dphidy*var_y + dphidz*dphidz*var_z;
        var_phi *= pow(180.0/M_PI,2.0); // convert to deg²
    }
}

