/***********************************************
 * Analysis of tracks	
 *
 * @author Felix Touchte Codjo
 * @date May 21, 2025
 * ********************************************/

/** Notes on the cuts
 * 0) requires at least one track, startTime > 0, one electron 
 * Exclude 
 * 1) electron -> |vz| > 15cm
 * 2) track -> theta < 5 
 *      - no number of hits == 2
 *      - always nhits from 4
 *      - always path at 0
 *      - vz for track always above 15 cm
 * 3) track -> |vz| > 15 cm
 *      - always nhits from 4
 *      - almost no longer path at 0
 * 4) track -> nhits >= 6
 * 5) Q2 < 14 and Q2 > 13.7
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>

#include "reader.h"
#include "futils.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TStyle.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

const double Ee = 2.23951; // incident energy if the electron, GeV
const double me = 0.511e-3; // energy mass of electron, GeV
const double Mp = 938.272e-3; // energy mass of proton, GeV
const double M_He = 3.73; // energy mass of Helium-4, GeV
const double Mt = M_He; // target rest mass : choose Mp or M_He, GeV

int main(int argc, char const *argv[]) {
		const char * filename = "/home/touchte-codjo/Desktop/hipofiles/occupancy/rec_clas_021317.evio.00000.hipo";
        hipo::reader  reader(filename);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
		hipo::bank  track0Bank(factory.getSchema("AHDC::track")); // to match old cooked files where dEdx was filled in AHDC::track instead of AHDC::kftrack
		hipo::bank  particleBank(factory.getSchema("REC::Particle"));
		hipo::bank  recEventBank(factory.getSchema("REC::Event"));
		hipo::event event;
		long unsigned int nevents = 0;
		long unsigned int nelectrons = 0;
		long unsigned int ntracks = 0;
		long unsigned int nentries = 0;
		// Histograms
		// electron
		TH1D* H1_vz_el    = new TH1D("vz_electron", "vz electron (cm)", 100, -16, 16); 
		TH1D* H1_p_el     = new TH1D("p_electron", "p electron (GeV)", 100, 2.17, 2.24); // 0.8, 2.17
		TH1D* H1_pT_el    = new TH1D("pT_electron", "pT_electron (GeV)", 100, 0.2, 0.5); // (0.07, 0.8)
		TH1D* H1_theta_el = new TH1D("theta_electron", "theta electron (deg)", 100, 5, 13); // (4, 25)
		TH1D* H1_phi_el   = new TH1D("phi_electron", "phi electron (deg)", 100, 0, 361);
		TH1D* H1_nelectrons_per_evt = new TH1D("nelectrons_per_evt", "# e^{-} / evt", 100, 0, 2);
		// ahdc (kf) track
		TH1D* H1_vz      = new TH1D("z_track", "z (cm)", 100, -16, 16); 
		TH1D* H1_p       = new TH1D("p_track", "p (GeV)", 100, 0.3, 3);
		TH1D* H1_pT      = new TH1D("pT_track", "pT_track (GeV)", 100, 0, 1.5);
		TH1D* H1_theta   = new TH1D("theta_track", "theta (deg)", 100, 0, 181);
		TH1D* H1_phi     = new TH1D("phi_track", "#phi (deg)", 100, 0, 361);
		TH1D* H1_nhits   = new TH1D("nhits", "number of hits per track", 100, 0, 15);
		TH1D* H1_adc     = new TH1D("sum_adc", "#Sigma #it{adc}", 100, 0, 25000);
		TH1D* H1_path    = new TH1D("sum_path", "path (mm)", 100, 0, 400);
		TH1D* H1_dEdx    = new TH1D("sum_dEdx", "dEdx (adc/mm)", 100, 0, 200);
		TH1D* H1_sum_res = new TH1D("sum_res", "#Sigma #it{res} (mm)", 100, -40, 3);
		TH1D* H1_chi2    = new TH1D("chi2", "chi2 (mm^{2})", 100, 0, 150);
		TH1D* H1_p_drift = new TH1D("p_drift", "p_{drift} (GeV)", 100, 0.3, 3);
		// pid ? correlation study
		TH2D* H2_p_dEdx    = new TH2D("p_drift_dEdx", "p_drift vs dEdx", 100, 0, 3, 100, 0, 200); H2_p_dEdx->Draw("COLZ");
		TH2D* H2_vze_vz    = new TH2D("vze_vz", "vz_{e} vs vz", 100, -16, 16, 100, -16, 16); H2_vze_vz->Draw("COLZ");
		TH1D* H1_delta_vz  = new TH1D("delta_vz", "#Delta vz = vz_{e} - vz", 100, -30, 30);
		TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta #phi = #phi_{e} - #phi", 100, -7, 7);
		// physics
		TH1D* H1_Q2 = new TH1D("Q2", "Q^{2} (GeV^{2})", 100, 0, 0.4);
		TH1D* H1_W2 = new TH1D("W2", "W^{2} (GeV^{2})", 100, 13.7, 14); // 24.8, 18, 14
		TH1D* H1_xB = new TH1D("xB", "x_{B}", 100, 0, 10); // (0, 1.2)
		TH1D* H1_nu = new TH1D("nu", "#nu = #Delta E = E - E' (GeV)", 100, 0, 0.06); // 1.42, 0.06

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
            // ----- Ignore event if
            // - no track
            // - startTime < 0 
            // - no electrons
            if ((trackBank.getRows() < 1) || 
                (recEventBank.get("startTime", 0) < 0)) { 
                continue;
            }
            int nelectrons_per_evt = 0;
			for (int i = 0; i < particleBank.getRows(); i++) {
				if (particleBank.getInt("pid", i) == 11) { // cut on electrons
                    nelectrons_per_evt++;
                }
            }
            if (nelectrons_per_evt == 0) continue;
            nelectrons += nelectrons_per_evt;
            H1_nelectrons_per_evt->Fill(nelectrons_per_evt);
            // ----- end Ingore  event
            // actually, we can only detect one electron!
            // physics and control cuts
            double Q2 = -9999;
            double nu = -9999;
            double W2 = -9999;
            double xB = -9999;
            double vz_e = -9999;
            double px_e = -999, py_e = -9999, pz_e = -9999;
            double p_e = -9999, theta_e = -9999, phi_e = -9999;
            // end physics
			for (int i = 0; i < particleBank.getRows(); i++) {
				if (particleBank.getInt("pid", i) == 11) { // electron
                    px_e = particleBank.getFloat("px",i);
                    py_e = particleBank.getFloat("py",i);
                    pz_e = particleBank.getFloat("pz",i);
					futils::cart2polar(px_e, py_e, pz_e, p_e, theta_e, phi_e);
                    vz_e = particleBank.getFloat("vz", i);
					Q2 = 4*Ee*sqrt(p_e*p_e + me*me)*pow(sin(theta_e/2),2); // Ee' ~ p as me << 1 GeV
					nu = Ee - p_e;
					W2 = Mt*Mt + 2*Mt*nu - Q2;
					xB = Q2/(2*Mp*nu);
				}
			}
            // CUT CAN BE DONE HERE
            if ((fabs(vz_e) > 15) ||
                (W2 > 14) || (W2 < 13.7)
            ) { continue;}
            nentries++; // number of electrons that pass the cut
            // Fill histograms for electrons
            H1_p_el->Fill(p_e); // here is already in GeV
            H1_pT_el->Fill(sqrt(pow(px_e,2) + pow(py_e,2))); // here is already in GeV
            H1_theta_el->Fill(theta_e*180.0/M_PI);
            H1_phi_el->Fill(phi_e*180/M_PI);
            H1_vz_el->Fill(vz_e);
            H1_Q2->Fill(Q2);
            H1_W2->Fill(W2);
            H1_nu->Fill(nu);
            H1_xB->Fill(xB);
			// track
			double  vz_first_track = -999;
            double phi_first_track = -999;
            int one = 0;
			for (int i = 0; i < trackBank.getRows(); i++) {
				double p, theta, phi;
				futils::cart2polar(trackBank.getFloat("px",i), trackBank.getFloat("py",i), trackBank.getFloat("pz",i), p, theta, phi);
                // CUT CAN BE DONE HERE
                double theta_deg = theta*180.0/M_PI;
                if (std::isnan(theta) || (theta_deg < 5) || 
                    (fabs(0.1*trackBank.getFloat("z", i)) > 15) ||
                    (trackBank.getInt("n_hits", i) < 6)
                ){ continue;}
				ntracks++;
                one++;
				if (one == 1) {
                    vz_first_track  = 0.1*trackBank.getFloat("z",i); // convert mm to cm
                    phi_first_track = phi;
                }
                H1_p->Fill(0.001*p); // here p was in MeV
				H1_pT->Fill(0.001*sqrt(pow(trackBank.getFloat("px",i),2) + pow(trackBank.getFloat("py",i),2))); 
				H1_theta->Fill(theta*180.0/M_PI);
				H1_phi->Fill(phi*180/M_PI);
				H1_vz->Fill(0.1*trackBank.getFloat("z", i)); // convert mm to cm
				H1_nhits->Fill(trackBank.getInt("n_hits", i)); 
				H1_adc->Fill(trackBank.getInt("sum_adc", i));
				H1_path->Fill(trackBank.getFloat("path", i));   
				H1_dEdx->Fill(track0Bank.getFloat("dEdx", i));    
				H1_sum_res->Fill(trackBank.getFloat("sum_residuals", i));
				H1_chi2->Fill(trackBank.getFloat("chi2", i));
				H1_p_drift->Fill(0.001*trackBank.getFloat("p_drift", i));
				H2_p_dEdx->Fill(0.001*trackBank.getFloat("p_drift", i), track0Bank.getFloat("dEdx", i));
			}
            if (one >= 1) { // require at least one track
                H2_vze_vz->Fill(vz_e, vz_first_track); 
                H1_delta_vz->Fill(vz_e - vz_first_track); 
                H1_delta_phi->Fill(phi_e - phi_first_track); 
            }
		}
		printf("nevents    : %ld \n", nevents);
		printf("nelectrons : %ld \n", nelectrons);
		printf("ntracks    : %ld \n", ntracks);
		printf("nentries   : %ld \n", nentries);
		// output
		const char * output = "./elastics.root";
		TFile *f = new TFile(output, "RECREATE");
        TDirectory *physics_dir  = f->mkdir("physics");
        TDirectory *electron_dir = f->mkdir("electron");
        TDirectory *track_dir    = f->mkdir("track");
        TDirectory *corr_dir     = f->mkdir("correlation");
            // physics
        physics_dir->cd();
		H1_Q2->Write("Q2");
		H1_W2->Write("W2");
		H1_xB->Write("xB");
		H1_nu->Write("nu");
            // electrons
        electron_dir->cd();
		H1_p_el->Write("p_el");
		H1_pT_el->Write("pT_el");
		H1_theta_el->Write("theta_el");
		H1_phi_el->Write("phi_el");
		H1_vz_el->Write("vz_el");
        H1_nelectrons_per_evt->Write("nelectrons_per_evt");
            // ahdc kf track
        track_dir->cd();
		H1_p->Write("p");
		H1_pT->Write("pT");
		H1_theta->Write("theta");
		H1_phi->Write("phi");
		H1_vz->Write("vz");
		H1_nhits->Write("nhits");
		H1_adc->Write("sum_adc");
		H1_path->Write("path");
		H1_dEdx->Write("dEdx");
		H1_sum_res->Write("sum_residuals");
		H1_chi2->Write("chi2");
		H1_p_drift->Write("p_drift");
            // correlations
        corr_dir->cd();
		H2_p_dEdx->Write("p_drift_dEdx");
		H2_vze_vz->Write("vz_electron_track");
		H1_delta_vz->Write("delta_vz");
		H1_delta_phi->Write("delta_phi");
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
