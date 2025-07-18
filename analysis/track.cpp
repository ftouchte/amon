/***********************************************
 * Analysis of tracks	
 *
 * @author Felix Touchte Codjo
 * @date May 21, 2025
 * ********************************************/

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
#include "TGraph.h"
#include "TStyle.h"

const double Ee = 10.6; // incident energy if the electron, GeV
const double me = 0.511e-3; // energy mass of electron, GeV
const double Mp = 938.272e-3; // energy mass of proton, GeV
const double M_He = 3.73; // energy mass of Helium-4, GeV
const double Mt = M_He; // target rest mass : choose Mp or M_He, GeV

int main(int argc, char const *argv[]) {
	if (argc > 1){
		hipo::reader  reader(argv[1]);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  track(factory.getSchema("AHDC::kftrack"));
		hipo::bank  track0(factory.getSchema("AHDC::track"));
		hipo::bank  particles(factory.getSchema("REC::Particle"));
		hipo::event event;
		long unsigned int nevents =0;
		long unsigned int nelectrons =0;
		long unsigned int ntracks =0;
		// Histograms
		// electron
		TH1D* H1_vz_el    = new TH1D("vz_electron", "vz electron (cm)", 100, -40, 40); 
		TH1D* H1_p_el     = new TH1D("p_electron", "p electron (GeV)", 100, 0, 11);
		TH1D* H1_theta_el = new TH1D("theta_electron", "theta electron (deg)", 100, 0, 180);
		TH1D* H1_phi_el   = new TH1D("phi_electron", "phi electron (deg)", 100, 0, 360);
		// ahdc (kf) track
		TH1D* H1_vz      = new TH1D("z_track", "z (cm)", 100, -40, 40); 
		TH1D* H1_p       = new TH1D("p_track", "p (GeV)", 100, 0, 11);
		TH1D* H1_theta   = new TH1D("theta_track", "theta (deg)", 100, 0, 180);
		TH1D* H1_phi     = new TH1D("phi_track", "phi (deg)", 100, 0, 360);
		TH1D* H1_nhits   = new TH1D("nhits", "number of hits per track", 100, 0, 15);
		TH1D* H1_adc     = new TH1D("sum_adc", "#Sigma #it{adc}", 100, 0, 40000);
		TH1D* H1_path    = new TH1D("sum_path", "path (mm)", 100, 0, 1000);
		TH1D* H1_dEdx    = new TH1D("sum_dEdx", "dEdx (adc/mm)", 100, 0, 500);
		TH1D* H1_sum_res = new TH1D("sum_res", "#Sigma #it{res} (mm)", 100, 0, 3);
		TH1D* H1_chi2    = new TH1D("chi2", "#Sigma #it{res} (mm^{2})", 100, 0, 1000);
		TH1D* H1_p_drift = new TH1D("p_drift", "p_{drift} (MeV ? in bankdef)", 100, 0, 11000);
		// pid ? correlation study
		TH2D* H2_dEdx_p  = new TH2D("dEdx_p_drift", "dEdx vs p_drift", 100, 0, 500, 100, 0, 11000); H2_dEdx_p->Draw("COLZ");
		TH2D* H2_vze_vz   = new TH2D("vze_vz", "vz_{e} vs vz", 100, -40, 40, 100, -40, 40); H2_vze_vz->Draw("COLZ");
		// physics
		TH1D* H1_Q2 = new TH1D("Q2", "Q^{2} (GeV^{2})", 100, 0, 40);
		TH1D* H1_W2 = new TH1D("W2", "W^{2} (GeV^{2})", 100, -20, 100);
		TH1D* H1_nu = new TH1D("nu", "#nu = #Delta E = E - E' (GeV)", 100, 0, 10.7);
		TH1D* H1_xB = new TH1D("xB", "x_{B}", 100, 0, 2);

		// Loop over events
		while( reader.next()){
			nevents++;
			//if (nevents > 10000) break;
			if (nevents % 100000 == 0) { printf("Start event %ld\n", nevents);}
			reader.read(event);
			event.getStructure(track);
			event.getStructure(track0);
			event.getStructure(particles);
			double vze = -999;
			// electrons
			for (int i = 0; i < particles.getRows(); i++) {
				int nelectrons_per_evt = 0;
				if (particles.getInt("pid", i) == 11) { // cut on electrons
					nelectrons++;
					nelectrons_per_evt++;
					double p, theta, phi;
					futils::cart2polar(particles.getFloat("px",i), particles.getFloat("py",i), particles.getFloat("pz",i), p, theta, phi);
					double Q2 = 4*Ee*sqrt(p*p + me*me)*sin(theta/2); // Ee' ~ p as me << 1 GeV
					double nu = Ee - p;
					double W2 = Mt*Mt + 2*Mt*nu - Q2;
					double xB = Q2/(2*Mt*nu);
					//printf("%lf ", xB);
					// elastic cuts : W2 ~ Mt
					//if ((W2 > 0.97) && (W2 < 1.03)) {
						H1_p_el->Fill(p); // here is already in GeV
						H1_theta_el->Fill(theta*180.0/M_PI);
						H1_phi_el->Fill(phi*180/M_PI);
						H1_vz_el->Fill(particles.getFloat("vz", i));
						H1_Q2->Fill(Q2);
						H1_W2->Fill(W2);
						H1_nu->Fill(nu);
						H1_xB->Fill(xB);
					//}
					if (nelectrons_per_evt == 1) { vze = particles.getFloat("vz", i);}
				}
			}
			// track
			for (int i = 0; i < track.getRows(); i++) {
				ntracks++;
				double p, theta, phi;
				futils::cart2polar(track.getFloat("px",i), track.getFloat("py",i), track.getFloat("pz",i), p, theta, phi);
				H1_p->Fill(p/1000.0); // here p was in MeV
				H1_theta->Fill(theta*180.0/M_PI);
				H1_phi->Fill(phi*180/M_PI);
				H1_vz->Fill(track.getFloat("z", i)/10.0); // convert mm to cm
				H1_nhits->Fill(track.getInt("n_hits", i)); 
				H1_adc->Fill(track.getInt("sum_adc", i));
				H1_path->Fill(track.getFloat("path", i));   
				H1_dEdx->Fill(track0.getFloat("dEdx", i));    
				H1_sum_res->Fill(track.getFloat("sum_residuals", i));
				H1_chi2->Fill(track.getFloat("chi2", i));
				H1_p_drift->Fill(track.getFloat("p_drift", i));
				H2_dEdx_p->Fill(track0.getFloat("dEdx", i), track.getFloat("p_drift", i));
			}
			H2_vze_vz->Fill(vze, track.getFloat("z", 0)/10.0); // first electrons vs the first track (en cm) 
		}
		printf("nevents    : %ld \n", nevents);
		printf("nelectrons : %ld \n", nelectrons);
		printf("ntracks    : %ld \n", ntracks);
		// output
		const char * output = "./output/track_study.root";
		TFile *f = new TFile(output, "RECREATE");
		H1_Q2->Write("Q2");
		H1_W2->Write("W2");
		H1_xB->Write("xB");
		H1_nu->Write("nu");
		H1_p_el->Write("p_el");
		H1_theta_el->Write("theta_el");
		H1_phi_el->Write("phi_el");
		H1_vz_el->Write("vz_el");
		H1_p->Write("p");
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
		H2_dEdx_p->Write("corr_dEdx_p_drift");
		H2_vze_vz->Write("corr_vze_vz");
		f->Close();
		printf("File created : %s\n", output);
	} 
	else {
		std::cout << " ***** please provide a filename..." << std::endl;
	}

	return 0;
}


