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

const double Ee = 10.6; // incident energy if the electron
const double me = 0.511e-3; // energy mass of electron
const double Mp = 938.272e-3; // energy mass of proton
const double M_He = 3.73; // energy mass of Helium-4
const double Mt = Mp; // target rest mass : choose Mp or M_He

int main(int argc, char const *argv[]) {
	if (argc > 1){
		hipo::reader  reader(argv[1]);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  track(factory.getSchema("AHDC::track"));
		hipo::bank  kftrack(factory.getSchema("AHDC::track"));
		hipo::bank  particles(factory.getSchema("REC::Particle"));
		hipo::event event;
		long unsigned int nevents =0;
		long unsigned int nelectrons =0;
		long unsigned int nkftracks =0;
		long unsigned int ntracks =0;
		// Histograms
		// electron
		TH1D* H1_vz_el    = new TH1D("vz_electron", "vz electron (cm)", 100, -40, 40); 
		TH1D* H1_p_el     = new TH1D("p_electron", "p electron (GeV)", 100, 0, 11);
		TH1D* H1_theta_el = new TH1D("theta_electron", "theta electron (deg)", 100, 0, 180);
		TH1D* H1_phi_el   = new TH1D("phi_electron", "phi electron (deg)", 100, 0, 360);
		// ahdc track
		TH1D* H1_vz     = new TH1D("z_track", "z (mm)", 100, -400, 400); 
		TH1D* H1_p     = new TH1D("p_track", "p (GeV)", 100, 0, 11);
		TH1D* H1_theta = new TH1D("theta_track", "theta (deg)", 100, 0, 180);
		TH1D* H1_phi   = new TH1D("phi_track", "phi (deg)", 100, 0, 360);
		TH1D* H1_nhits = new TH1D("nhits", "number of hits per track", 100, 0, 15);
		TH1D* H1_adc   = new TH1D("sum_adc", "sum of adc", 100, 0, 40000);
		// ahdc kftrack
		TH1D* H1_vz_kf     = new TH1D("z_kftrack", "z (mm)", 100, -400, 400); 
		TH1D* H1_p_kf     = new TH1D("p_kftrack", "p (GeV)", 100, 0, 11);
		TH1D* H1_theta_kf = new TH1D("theta_kftrack", "theta (deg)", 100, 0, 180);
		TH1D* H1_phi_kf   = new TH1D("phi_kftrack", "phi (deg)", 100, 0, 360);
		// physics
		TH1D* H1_Q2 = new TH1D("Q2", "Q^{2} (GeV^{2})", 100, -100, 100);
		TH1D* H1_W2 = new TH1D("W2", "W^{2} (GeV^{2})", 100, -100, 100);
		TH1D* H1_nu = new TH1D("nu", "#nu = #Delta E = E - E' (GeV)", 100, 0, 11);
		TH1D* H1_xB = new TH1D("xB", "x_{B}", 100, 0, 1.05);

		// Loop over events
		while( reader.next()){
			nevents++;
			//if (nevents > 10000) break;
			if (nevents % 100000 == 0) { printf("Start event %ld\n", nevents);}
			reader.read(event);
			event.getStructure(track);
			event.getStructure(kftrack);
			event.getStructure(particles);
			// electrons
			for (int i = 0; i < particles.getRows(); i++) {
				if (particles.getInt("pid", i) == 11) { // cut on electrons
					nelectrons++;
					double p, theta, phi;
					futils::cart2polar(particles.getFloat("px",i), particles.getFloat("py",i), particles.getFloat("pz",i), p, theta, phi);
					double Q2 = 4*Ee*p*sin(theta/2); // Ee' ~ p
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
				H1_vz->Fill(track.getFloat("z", i));
			}
			// kftrack
			for (int i = 0; i < kftrack.getRows(); i++) {
				nkftracks++;
				double p, theta, phi;
				futils::cart2polar(kftrack.getFloat("px",i), kftrack.getFloat("py",i), kftrack.getFloat("pz",i), p, theta, phi);
				H1_p_kf->Fill(p/1000.0); // here p was in MeV
				H1_theta_kf->Fill(theta*180.0/M_PI);
				H1_phi_kf->Fill(phi*180/M_PI);
				H1_vz_kf->Fill(kftrack.getFloat("z", i));
				H1_nhits->Fill(kftrack.getInt("n_hits", i));
				H1_adc->Fill(kftrack.getInt("sum_adc", i));
			}
		}
		printf("nevents    : %ld \n", nevents);
		printf("nelectrons : %ld \n", nelectrons);
		printf("nkftracks  : %ld \n", nkftracks);
		printf("ntracks    : %ld \n", nkftracks);
		// output
		const char * output = "../output/proton_track_study.root";
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
		H1_p_kf->Write("p_kf");
		H1_theta_kf->Write("theta_kf");
		H1_phi_kf->Write("phi_kf");
		H1_vz_kf->Write("vz_kf");
		f->Close();
		printf("File created : %s\n", output);
	} 
	else {
		std::cout << " ***** please provide a filename..." << std::endl;
	}

	return 0;
}


