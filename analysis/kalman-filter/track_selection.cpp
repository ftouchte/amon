/***************************************************
 * Small code to determine the cuts to be applied
 * to select good track
 * 
 * @author ftouchte
 * @date April 21, 2026
 ***************************************************/

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
#include "TLine.h"

#include "futils.h"
#include "fOptions.h"
#include "Units.h"

void progressBar(int state, int bar_length = 100);
TCanvas* showCuts(TH1D* h, double lim_inf, double lim_sup);

/// --- Constants
const double beam_energy = 2.23951 * Units::GeV; // incident energy if the electron, GeV
const double electron_mass = 0.511e-3 * Units::GeV; // energy mass of electron, GeV
const double proton_mass = 938.272e-3 * Units::GeV; // energy mass of proton, GeV
const double helium_mass = 3.73 * Units::GeV; // energy mass of Helium-4, GeV
const double deuteron_mass = 1.875 * Units::GeV; // energy mass of Deuterium, GeV

int main(int argc, char const *argv[]) {

    /// --- Selection cuts
    double chi2_min = 0.1;
    double chi2_max = 3;
    double vz_min = -24; //cm
    double vz_max = 15; // cm
    double nhits_min = 7;

    /// --- start timer
    auto start = std::chrono::high_resolution_clock::now();

    /// --- Load options
    fOptions OPT({"-i", "-o", "-v", "-simu"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    /// --- Initialise histograms
    TH2D* H2_corr_p_dEdx = new TH2D("corr_p_dEdx", "dEdx vs p; p (MeV); dEdx (MeV/mm)", 100, 80, 600, 100, 0, 220);
    TH2D* H2_corr_p_dEdx_selection = new TH2D("corr_p_dEdx_selection", "dEdx vs p; p (MeV); dEdx (MeV/mm)", 100, 80, 600, 100, 0, 220);
    TH1D* H1_track_vz = new TH1D("track_vz", "vz; vz (cm); count", 100, -30, 25);
    TH1D* H1_track_nhits = new TH1D("track_nhits", "nhits; nhits; count", 20, 0, 20);
    TH1D* H1_track_p = new TH1D("track_p", "track p; p (MeV); count", 100, 80, 600);
    TH1D* H1_track_p_selection = new TH1D("track_p_selection", "track p; p (MeV); count", 100, 80, 600);
    TH1D* H1_track_theta = new TH1D("track_theta", "track theta; theta (deg); count", 100, 0, 180);
    TH1D* H1_track_theta_selection = new TH1D("track_theta_selection", "track theta; theta (deg); count", 100, 0, 180);
    TH1D* H1_track_phi = new TH1D("track_phi", "track phi; phi (deg); count", 100, 0, 360);
    TH1D* H1_track_phi_selection = new TH1D("track_phi_selection", "track phi; phi (deg); count", 100, 0, 360);
    TH1D* H1_track_chi2 = new TH1D("track_chi2", "track chi2; chi2; count", 100, 0, 8);
    TH1D* H1_W2 = new TH1D("W2", "W^{2}; W^{2} (GeV^{2}); count", 100, 3.2, 6);

    /// --- Start analysis
    int nfile = 0;
    long unsigned int nevents =0;
    long unsigned int ntracks =0;
    for (auto file : filenames) {
        nfile++;
        printf("> Open file %d/%d : %s\n", nfile, (int) filenames.size(), file.c_str());
        
        /// --- Initialise HIPO reader
        hipo::reader  reader(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        /// --- Defien banks to be read
        hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
        hipo::bank  recBank(factory.getSchema("REC::Particle"));


        /// --- Loop over events
        hipo::event event;
        long unsigned int nevents_per_file = 0;

        while( reader.next()){
            nevents++;
            nevents_per_file++;

            // display progress Bar
            if ((nevents_per_file % 1000 == 0) || ((int) nevents_per_file == reader.getEntries())) {
                progressBar(100.0*nevents_per_file/reader.getEntries());
            }
            // load bank content for this event
            reader.read(event);
            event.getStructure(trackBank);
            event.getStructure(recBank);

            bool aGoodTrackHasBeennFound = false;
            for (int row = 0; row < trackBank.getRows(); row++) 
            {
                double px = trackBank.getFloat("px", row); // MeV
                double py = trackBank.getFloat("py", row); // MeV
                double pz = trackBank.getFloat("pz", row); // MeV
                double p = sqrt(pow(px, 2)+pow(py, 2)+pow(pz, 2)); // MeV
                double dEdx = trackBank.getFloat("dEdx", row);
                int nhits = trackBank.getInt("n_hits", row);
                double vz = trackBank.getFloat("z", row)*0.1; // cm
                double chi2 = trackBank.getFloat("chi2", row);

                double theta_deg = acos(pz/p)*180.0/M_PI;
                double phi_rad = atan2(py, px);
                if (phi_rad < 0) phi_rad += 2*M_PI;
                double phi_deg = phi_rad*180.0/M_PI;

                //boolean flag = nhits >= 7 && p > 100 && p < 600 && dEdx > 0 && dEdx < 200;
                // if (nhits >= 7 && p > 100 && p < 600 && dEdx > 0 && dEdx < 200) {
                //     trackRows.add(row);
                // }

                H2_corr_p_dEdx->Fill(p, dEdx);
                H1_track_vz->Fill(vz);
                H1_track_nhits->Fill(nhits);
                H1_track_p->Fill(p);
                H1_track_phi->Fill(phi_deg);
                H1_track_theta->Fill(theta_deg);
                H1_track_chi2->Fill(chi2);

                if (nhits >= nhits_min && vz > vz_min && vz < vz_max && chi2 > chi2_min && chi2 < chi2_max) {
                    H2_corr_p_dEdx_selection->Fill(p, dEdx);
                    H1_track_p_selection->Fill(p);
                    H1_track_phi_selection->Fill(phi_deg);
                    H1_track_theta_selection->Fill(theta_deg);
                    ntracks++;
                    aGoodTrackHasBeennFound = true;
                }

            } // end loop over track rows

            if (aGoodTrackHasBeennFound) {
                for (int row = 0; row < recBank.getRows(); row++) {
                    // Select trigger electrons
                    if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0) {
                        // compute kinematic variables
                        double px = recBank.getFloat("px", row) * Units::GeV;
                        double py = recBank.getFloat("py", row) * Units::GeV;
                        double pz = recBank.getFloat("pz", row) * Units::GeV;
                        double p = sqrt(pow(px, 2)+pow(py, 2)+pow(pz, 2));
                        // physics kinematics
                        double scattered_beam_energy = sqrt(pow(p,2) + pow(electron_mass,2));
                        double nu = beam_energy - scattered_beam_energy;
                        double theta = acos(pz/p) * Units::rad;
                        double Q2 = 4*beam_energy*scattered_beam_energy*pow(sin(theta/2),2);
                        double W2 = pow(deuteron_mass,2) + 2*deuteron_mass*nu - Q2;
                        H1_W2->Fill(W2);

                    } // end selection of trigger electron
                } // end loop over recBank rows
            } // end condition on good track

        } // end loop over events for this file
        printf("\033[34m nevents : %ld \033[0m \n", nevents);


    } // end loop over input files

    

    /// Store histograms
    TFile *f = new TFile(output.c_str(), "RECREATE");
    
    H1_track_chi2->Write("track_chi2");
    showCuts(H1_track_chi2, chi2_min, chi2_max)->Write("track_chi2_with_cuts");

    H1_track_vz->Write("track_vz");
    showCuts(H1_track_vz, vz_min, vz_max)->Write("track_vz_with_cuts");

    H1_track_nhits->Write("track_nhits");
    showCuts(H1_track_nhits, nhits_min, 1e10)->Write("track_nhits_with_cuts");

    H2_corr_p_dEdx->Write("corr_p_dEdx");
    H2_corr_p_dEdx_selection->Write("corr_p_dEdx_with_cuts");

    H1_track_p->Write("track_p");
    H1_track_p_selection->Write("track_p_with_cuts");

    H1_track_theta->Write("track_theta");
    H1_track_theta_selection->Write("track_theta_with_cuts");

    H1_track_phi->Write("track_phi");
    H1_track_phi_selection->Write("track_phi_with_cuts");

    H1_W2->Write("W2");
    //showCuts(H1_W2, 3.5, 3.8)->Write("W2_showing_elastic_limits");

    f->Close();
    printf("nb good tracks : %ld\n", ntracks);
    printf("File created : %s\n", output.c_str());

    /// --- end of the program
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

TCanvas* showCuts(TH1D* h, double lim_inf, double lim_sup) {
    TCanvas* canvas = new TCanvas();
    canvas->cd();
    h->Draw();
    double ymax = h->GetMaximum();
    TLine* line_inf = new TLine(lim_inf, 0, lim_inf, ymax);
    TLine* line_sup = new TLine(lim_sup, 0, lim_sup, ymax);
    
    line_inf->SetLineColor(kRed);
    line_sup->SetLineColor(kRed);
    line_inf->SetLineWidth(2);
    line_sup->SetLineWidth(2);
    
    line_inf->Draw("same L");
    line_sup->Draw("same L");

    return canvas;
}