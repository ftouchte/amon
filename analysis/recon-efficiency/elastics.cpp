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
 * 5) Q2 < 14 and Q2 > 13.85
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
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

const double Ee = 2.23951; // incident energy if the electron, GeV
const double me = 0.511e-3; // energy mass of electron, GeV
const double Mp = 938.272e-3; // energy mass of proton, GeV
const double M_He = 3.73; // energy mass of Helium-4, GeV
const double Mt = M_He; // target rest mass : choose Mp or M_He, GeV

int main(int argc, char const *argv[]) {
		//const char * filename = "/home/touchte-codjo/Desktop/hipofiles/occupancy/rec_clas_021317.evio.00000.hipo";
		const char * filename = "/home/touchte-codjo/Desktop/hipofiles/track/He4/all_rec_clas_021317.hipo";
        double W2_min = 13.85;
        double W2_max = 14;
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
        TCanvas* canvas = new TCanvas("Canvas_1");
		// electron
		TH1D* H1_vz_el    = new TH1D("vz_electron", "vz electron (cm)", 100, -40, 40); 
		TH1D* H1_vz_el_all    = new TH1D("vz_electron_all", "vz electron (cm)", 100, -40, 40); 
		TH1D* H1_p_el     = new TH1D("p_electron", "p electron (GeV)", 100, 2.17, 2.24); // 0.8, 2.17
		TH1D* H1_pT_el    = new TH1D("pT_electron", "pT_electron (GeV)", 100, 0.2, 0.45); // (0.07, 0.8)
		TH1D* H1_theta_el = new TH1D("theta_electron", "#theta electron (deg)", 100, 4, 25); // (4, 25)
		TH1D* H1_phi_el   = new TH1D("phi_electron", "#phi electron (deg)", 100, 0, 361);
		TH1D* H1_nelectrons_per_evt = new TH1D("nelectrons_per_evt", "# e^{-} / evt", 100, 0, 2);
		// ahdc (kf) track
		TH1D* H1_vz      = new TH1D("z_track", "z (cm)", 100, -30, 30); 
		TH1D* H1_vz_all      = new TH1D("z_track_all", "z (cm)", 100, -30, 30); 
		TH1D* H1_p       = new TH1D("p_track", "p (GeV)", 100, 0.3, 3);
		TH1D* H1_pT      = new TH1D("pT_track", "pT_track (GeV)", 100, 0.2, 0.73);
		TH1D* H1_theta   = new TH1D("theta_track", "#theta (deg)", 100, 0, 181);
		TH1D* H1_theta_all   = new TH1D("theta_track_all", "#theta (deg)", 100, 0, 181);
		TH1D* H1_phi     = new TH1D("phi_track", "#phi (deg)", 100, 0, 361);
		TH1D* H1_nhits   = new TH1D("nhits", "number of hits per track", 100, 0, 15);
		TH1D* H1_nhits_all   = new TH1D("nhits_all", "number of hits per track", 100, 0, 15);
		TH1D* H1_adc     = new TH1D("sum_adc", "#Sigma #it{adc}", 100, 0, 25000);
		TH1D* H1_path    = new TH1D("sum_path", "path (mm)", 100, 40, 300);
		TH1D* H1_dEdx    = new TH1D("sum_dEdx", "dEdx (adc/mm)", 100, 0, 400);
		TH1D* H1_sum_res = new TH1D("sum_res", "#Sigma #it{res} (mm)", 100, -25, 3);
		TH1D* H1_sum_res_ndef = new TH1D("sum_res_ndef", "#Sigma #it{res} / nhits (mm)", 100, -4.17, 0.5);
		TH1D* H1_chi2    = new TH1D("chi2", "chi2 (mm^{2})", 100, 0, 100);
		TH1D* H1_chi2ndef    = new TH1D("chi2ndef", "chi2/ndef (mm^{2})", 100, 0, 16.67);
		TH1D* H1_p_drift = new TH1D("p_drift", "p_{drift} (GeV)", 100, 0.3, 3);
		// pid ? correlation study
		TH2D* H2_p_dEdx    = new TH2D("pTe_dEdx", "pT electron vs dEdx", 100, 0.2, 0.45, 100, 0, 400); H2_p_dEdx->Draw("COLZ");
		TH2D* H2_vze_vz    = new TH2D("vze_vz", "vz_{e} vs vz", 100, -25, 10, 100, -16, 16); H2_vze_vz->Draw("COLZ");
		TH2D* H2_pTe_pT    = new TH2D("pTe_pT", "pT_{e} vs pT", 100, 0.2, 0.45, 100, 0.2, 0.73); H2_pTe_pT->Draw("COLZ");
		TH1D* H1_delta_vz  = new TH1D("delta_vz", "#Delta vz = vz_{e} - vz (cm)", 100, -30, 30);
		TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta #phi = #phi_{e} - #phi (deg)", 100, -260, 200);
		// physics
		TH1D* H1_Q2 = new TH1D("Q2", "Q^{2} (GeV^{2})", 100, 0, 0.4);
		TH1D* H1_W2 = new TH1D("W2", "W^{2} (GeV^{2})", 100, 13.7, 14); // 24.8, 18, 14
		TH1D* H1_W2_all = new TH1D("W2_all", "W^{2} (GeV^{2})", 100, 13.7, 18); // 24.8, 18, 14
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
            double pT_e = -9999;
            // end physics
			for (int i = 0; i < particleBank.getRows(); i++) {
				if (particleBank.getInt("pid", i) == 11) { // electron
                    px_e = particleBank.getFloat("px",i);
                    py_e = particleBank.getFloat("py",i);
                    pz_e = particleBank.getFloat("pz",i);
					futils::cart2polar(px_e, py_e, pz_e, p_e, theta_e, phi_e);
                    pT_e = sqrt(pow(px_e,2) + pow(py_e,2));
                    vz_e = particleBank.getFloat("vz", i);
					Q2 = 4*Ee*sqrt(p_e*p_e + me*me)*pow(sin(theta_e/2),2); // Ee' ~ p as me << 1 GeV
					nu = Ee - p_e;
					W2 = Mt*Mt + 2*Mt*nu - Q2;
					xB = Q2/(2*Mp*nu);
				}
			}
            H1_W2_all->Fill(W2);
            H1_vz_el_all->Fill(vz_e);
            // CUT CAN BE DONE HERE
            if ((vz_e < -25) || (vz_e > 10) ||
                (W2 > W2_max) || (W2 < W2_min)
            ) { continue;}
            nentries++; // number of electrons that pass the cut
            // Fill histograms for electrons
            H1_p_el->Fill(p_e); // here is already in GeV
            H1_pT_el->Fill(pT_e); // here is already in GeV
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
            double dEdx_first_track = -999;
            double pT_first_track = -999;
            int one = 0;
			for (int i = 0; i < trackBank.getRows(); i++) {
				double p, theta, phi;
				futils::cart2polar(trackBank.getFloat("px",i), trackBank.getFloat("py",i), trackBank.getFloat("pz",i), p, theta, phi);
				H1_theta_all->Fill(theta*180.0/M_PI);
				H1_vz_all->Fill(0.1*trackBank.getFloat("z", i)); // convert mm to cm
				H1_nhits_all->Fill(trackBank.getInt("n_hits", i)); 
                // CUT CAN BE DONE HERE
                double theta_deg = theta*180.0/M_PI;
                if (std::isnan(theta) || (theta_deg < 10) || (theta_deg > 170) || 
                    (fabs(0.1*trackBank.getFloat("z", i)) > 15) ||
                    (trackBank.getInt("n_hits", i) < 6)
                ){ continue;}
				ntracks++;
                one++;
				if (one == 1) {
                    vz_first_track  = 0.1*trackBank.getFloat("z",i); // convert mm to cm
                    phi_first_track = phi;
                    pT_first_track = 0.001*sqrt(pow(trackBank.getFloat("px",i),2) + pow(trackBank.getFloat("py",i),2));
                    dEdx_first_track = track0Bank.getFloat("dEdx", i);
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
				H1_sum_res_ndef->Fill(trackBank.getFloat("sum_residuals", i)/trackBank.getInt("n_hits", i));
				H1_chi2->Fill(trackBank.getFloat("chi2", i));
				H1_chi2ndef->Fill(trackBank.getFloat("chi2", i)/trackBank.getInt("n_hits", i));
				H1_p_drift->Fill(0.001*trackBank.getFloat("p_drift", i));
			}
            if (one >= 1) { // require at least one track
                H2_vze_vz->Fill(vz_e, vz_first_track); 
                H1_delta_vz->Fill(vz_e - vz_first_track); 
                H1_delta_phi->Fill((phi_e - phi_first_track)*180/M_PI); // convert rad in deg 
				H2_p_dEdx->Fill(pT_e, dEdx_first_track);
                H2_pTe_pT->Fill(pT_e, pT_first_track);
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
        // W2 with a display of the cuts
            canvas->Clear();
            H1_W2_all->Draw();
            int w2_max = H1_W2_all->GetBinContent(H1_W2_all->GetMaximumBin());
            TGraph* gr_w2_min = new TGraph();
            gr_w2_min->AddPoint(13.85, 0);
            gr_w2_min->AddPoint(13.85, w2_max);
            gr_w2_min->SetLineColor(kRed);
            gr_w2_min->Draw("L");
            TGraph* gr_w2_max = new TGraph();
            gr_w2_max->AddPoint(14, 0);
            gr_w2_max->AddPoint(14, w2_max);
            gr_w2_max->SetLineColor(kRed);
            gr_w2_max->Draw("L");
            TLegend* legend_w2 = new TLegend(0.1,0.7,0.48,0.9);
            legend_w2->SetHeader("Cuts","C");
            legend_w2->AddEntry(gr_w2_min, (std::string("W2 >= ") + std::to_string(W2_min)).c_str());
            legend_w2->AddEntry(gr_w2_max, (std::string("W2 >= ") + std::to_string(W2_max)).c_str());
            legend_w2->Draw();
            canvas->Write("W2_cut_displayed");
		H1_xB->Write("xB");
		H1_nu->Write("nu");
        // electrons
        electron_dir->cd();
		H1_p_el->Write("p_el");
		H1_pT_el->Write("pT_el");
		H1_theta_el->Write("theta_el");
		H1_phi_el->Write("phi_el");
		H1_vz_el->Write("vz_el");
            // cut vz electron
            canvas->Clear();
            H1_vz_el_all->Draw();
            int vz_el_max = H1_vz_el_all->GetBinContent(H1_vz_el_all->GetMaximumBin());
            TGraph* gr_vz_el_min = new TGraph();
            gr_vz_el_min->AddPoint(-25, 0);
            gr_vz_el_min->AddPoint(-25, vz_el_max);
            gr_vz_el_min->SetLineColor(kRed);
            gr_vz_el_min->Draw("L");
            TGraph* gr_vz_el_max = new TGraph();
            gr_vz_el_max->AddPoint(10, 0);
            gr_vz_el_max->AddPoint(10, vz_el_max);
            gr_vz_el_max->SetLineColor(kRed);
            gr_vz_el_max->Draw("L");
            TLegend* legend_vz_el = new TLegend(0.1,0.7,0.48,0.9);
            legend_vz_el->SetHeader("Cuts","C");
            legend_vz_el->AddEntry(gr_vz_el_min, "vz >= -25");
            legend_vz_el->AddEntry(gr_vz_el_max, "vz <= 10");
            legend_vz_el->Draw();
            canvas->Write("vz_el_cut_displayed");
        H1_nelectrons_per_evt->Write("nelectrons_per_evt");
        // ahdc kf track
        track_dir->cd();
		H1_p->Write("p");
		H1_pT->Write("pT");
		H1_theta->Write("theta");
            // cut theta 
            canvas->Clear();
            H1_theta_all->Draw();
            int theta_max = H1_theta_all->GetBinContent(H1_theta_all->GetMaximumBin());
            TGraph* gr_theta_min = new TGraph();
            gr_theta_min->AddPoint(10, 0);
            gr_theta_min->AddPoint(10, theta_max);
            gr_theta_min->SetLineColor(kRed);
            gr_theta_min->Draw("L");
            TGraph* gr_theta_max = new TGraph();
            gr_theta_max->AddPoint(170, 0);
            gr_theta_max->AddPoint(170, theta_max);
            gr_theta_max->SetLineColor(kRed);
            gr_theta_max->Draw("L");
            TLegend* legend_theta = new TLegend(0.1,0.7,0.48,0.9);
            legend_theta->SetHeader("Cuts","C");
            legend_theta->AddEntry(gr_theta_min, "theta >= 10 deg");
            legend_theta->AddEntry(gr_theta_max, "theta <= 170 deg");
            legend_theta->Draw();
            canvas->Write("theta_cut_displayed");
		H1_phi->Write("phi");
		H1_vz->Write("vz");
            // cut vz 
            canvas->Clear();
            H1_vz_all->Draw();
            int vz_max = H1_vz_all->GetBinContent(H1_vz_all->GetMaximumBin());
            TGraph* gr_vz_min = new TGraph();
            gr_vz_min->AddPoint(-15, 0);
            gr_vz_min->AddPoint(-15, vz_max);
            gr_vz_min->SetLineColor(kRed);
            gr_vz_min->Draw("L");
            TGraph* gr_vz_max = new TGraph();
            gr_vz_max->AddPoint(15, 0);
            gr_vz_max->AddPoint(15, vz_max);
            gr_vz_max->SetLineColor(kRed);
            gr_vz_max->Draw("L");
            TLegend* legend_vz = new TLegend(0.1,0.7,0.48,0.9);
            legend_vz->SetHeader("Cuts","C");
            legend_vz->AddEntry(gr_vz_min, "vz >= -15");
            legend_vz->AddEntry(gr_vz_max, "vz <= 15");
            legend_vz->Draw();
            canvas->Write("vz_cut_displayed");
		H1_nhits->Write("nhits");
            // cut nhits 
            canvas->Clear();
            H1_nhits_all->Draw();
            int nhits_max = H1_nhits_all->GetBinContent(H1_nhits_all->GetMaximumBin());
            TGraph* gr_nhits_min = new TGraph();
            gr_nhits_min->AddPoint(6, 0);
            gr_nhits_min->AddPoint(6, nhits_max);
            gr_nhits_min->SetLineColor(kRed);
            gr_nhits_min->Draw("L");
            TLegend* legend_nhits = new TLegend(0.1,0.7,0.48,0.9);
            legend_nhits->SetHeader("Cuts","C");
            legend_nhits->AddEntry(gr_nhits_min, "nhits >= 6");
            legend_nhits->Draw();
            canvas->Write("nhits_cut_displayed");
		H1_adc->Write("sum_adc");
		H1_path->Write("path");
		H1_dEdx->Write("dEdx");
		H1_sum_res->Write("sum_residuals");
		H1_sum_res_ndef->Write("sum_residuals_per_nhits");
		H1_chi2->Write("chi2");
		H1_chi2ndef->Write("chi2ndef");
		H1_p_drift->Write("p_drift");
        // correlations
        corr_dir->cd();
        H2_p_dEdx->GetXaxis()->SetTitle("pT electron (GeV)"); 
        H2_p_dEdx->GetYaxis()->SetTitle("dEdx (adc/mm)");
		H2_p_dEdx->Write("pTe_dEdx");
        H2_vze_vz->GetXaxis()->SetTitle("vz electron (mm)"); 
        H2_vze_vz->GetYaxis()->SetTitle("vz track (mm)");
		H2_vze_vz->Write("vz_electron_track"); 
		H1_delta_vz->Write("delta_vz");
		H1_delta_phi->Write("delta_phi");
        // pT_e and pT
        H2_pTe_pT->GetXaxis()->SetTitle("pT electron (GeV)"); 
        H2_pTe_pT->GetYaxis()->SetTitle("pT track(GeV)");
		H2_pTe_pT->Write("pTe_pT");
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
