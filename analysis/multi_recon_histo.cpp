/***********************************************
 * Recontruction analysis	
 *
 * @author Felix Touchte Codjo
 * @date April 07, 2025
 * ********************************************/

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>

#include "reader.h"
#include "futils.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"

int main(int argc, char const *argv[])
{
    if (argc > 1){
	// AHDC::Track
	double px1,py1,pz1;
	double p1, theta1, phi1;
	double vz1;
	TH1D* hist_p1 	  = new TH1D("hist_p1","AHDC::Track #rightarrow p",100, 0,1000); // Mev
	TH1D* hist_theta1 = new TH1D("hist_theta1","AHDC::Track #rightarrow #theta",100, 0,190); // deg
	TH1D* hist_phi1   = new TH1D("hist_phi1","AHDC::Track #rightarrow #phi",100, 0,370); // phi
	TH1D* hist_vz1    = new TH1D("hist_vz1","AHDC::Track #rightarrow vz",100, 0, 350); // mm
	// REC::Particle
	double px2,py2,pz2;
	double p2, theta2, phi2;
	double vz2;
	TH1D* hist_p2     = new TH1D("hist_p2","REC::Particle, pid = 11 #rightarrow p",100, 0,1000); // Mev
	TH1D* hist_theta2 = new TH1D("hist_theta2","REC::Particle, pid = 11 #rightarrow #theta",100, 0,190); // deg
	TH1D* hist_phi2   = new TH1D("hist_phi2","REC::Particle, pid = 11 #rightarrow #phi",100, 0,370); // phi
	TH1D* hist_vz2    = new TH1D("hist_vz2","REC::Particle, pid = 11 #rightarrow vz",100, 0, 350); // mm
	// Others
	int nelectrons = 0;
	long unsigned int nEvent =0;
	for (int num_file = 1; num_file <= argc; num_file++) {
		std::cout << " >>>>>>>>>>>     file       :  " << argv[num_file] << std::endl;
		hipo::reader  reader(argv[num_file]);
		hipo::banklist banklist = reader.getBanks({"AHDC::Track","AHDC::adc","ATOF::tdc", "REC::Particle"});

		while( reader.next(banklist)){
			// AHDC::Track
			for(int col = 0; col < banklist[0].getRows(); col++){ 
				px1 = banklist[0].getFloat("px",col);
				py1 = banklist[0].getFloat("py",col);
				pz1 = banklist[0].getFloat("pz",col);
				vz1 = banklist[0].getFloat("z",col);
				/*p*/
				futils::cart2polar(px1,py1,pz1,p1,theta1,phi1);
				hist_p1->Fill(p1*1.0); // MeV
				hist_theta1->Fill(theta1*180/M_PI); // deg
				hist_phi1->Fill(phi1*180/M_PI); // deg
				hist_vz1->Fill(vz1);
			}
			// REC::Parcticle
			for(int col = 0; col < banklist[3].getRows(); col++){
				if (banklist[3].getInt("pid",col) == 11) {
					nelectrons++;
					px2 = banklist[0].getFloat("px",col);
					py2 = banklist[0].getFloat("py",col);
					pz2 = banklist[0].getFloat("pz",col);
					vz2 = banklist[0].getFloat("z",col);
					futils::cart2polar(px2,py2,pz2,p2,theta2,phi2);
					hist_p2->Fill(p2*1000.0); // MeV
					hist_theta2->Fill(theta2*180/M_PI); // deg
					hist_phi2->Fill(phi2*180/M_PI); // deg
					hist_vz2->Fill(vz2);
				}
			}
			// END	
			nEvent++;
		}
	}
        std::cout << "nEvent : " << nEvent << std::endl;
        std::cout << "nelectrons : " << nelectrons << std::endl;

        // *****************************************
        //          PlOT AHDC::Track
        // *****************************************

        TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
        canvas1->Divide(2,2);
        gStyle->SetOptStat("nemruo"); 
        // hist_p
        canvas1->cd(1);
        hist_p1->GetXaxis()->SetTitle("p (MeV)");
        hist_p1->GetXaxis()->SetTitleSize(0.05);
        hist_p1->GetYaxis()->SetTitle("#events");
        hist_p1->GetYaxis()->SetTitleSize(0.05);
        hist_p1->Draw();
        //hist_theta
        canvas1->cd(2);
        hist_theta1->GetXaxis()->SetTitle("#theta (deg)");
        hist_theta1->GetXaxis()->SetTitleSize(0.05);
        hist_theta1->GetYaxis()->SetTitle("#events");
        hist_theta1->GetYaxis()->SetTitleSize(0.05);
        hist_theta1->Draw();
        //hist_phi
        canvas1->cd(3);
        hist_phi1->GetXaxis()->SetTitle("#phi (deg)");
        hist_phi1->GetXaxis()->SetTitleSize(0.05);
        hist_phi1->GetYaxis()->SetTitle("#events");
        hist_phi1->GetYaxis()->SetTitleSize(0.05);
        hist_phi1->Draw();
        // hist_vz
        canvas1->cd(4);
        hist_vz1->GetXaxis()->SetTitle("vz (mm)");
        hist_vz1->GetXaxis()->SetTitleSize(0.05);
        hist_vz1->GetYaxis()->SetTitle("#events");
        hist_vz1->GetYaxis()->SetTitleSize(0.05);
        hist_vz1->Draw();
        // SAVE
        canvas1->Print("../output/ahdc_track.pdf");
        delete hist_p1; delete hist_theta1; delete hist_phi1; delete hist_vz1;
        delete canvas1;

        // *****************************************
        //          PlOT REC::Particle
        // *****************************************

        TCanvas* canvas2 = new TCanvas("c2","c2 title",1366,768);
        canvas2->Divide(2,2);
        gStyle->SetOptStat("nemruo"); 
        // hist_p
        canvas2->cd(1);
        hist_p2->GetXaxis()->SetTitle("p (MeV)");
        hist_p2->GetXaxis()->SetTitleSize(0.05);
        hist_p2->GetYaxis()->SetTitle("#events");
        hist_p2->GetYaxis()->SetTitleSize(0.05);
        hist_p2->Draw();
        //hist_theta
        canvas2->cd(2);
        hist_theta2->GetXaxis()->SetTitle("#theta (deg)");
        hist_theta2->GetXaxis()->SetTitleSize(0.05);
        hist_theta2->GetYaxis()->SetTitle("#events");
        hist_theta2->GetYaxis()->SetTitleSize(0.05);
        hist_theta2->Draw();
        //hist_phi
        canvas2->cd(3);
        hist_phi2->GetXaxis()->SetTitle("#phi (deg)");
        hist_phi2->GetXaxis()->SetTitleSize(0.05);
        hist_phi2->GetYaxis()->SetTitle("#events");
        hist_phi2->GetYaxis()->SetTitleSize(0.05);
        hist_phi2->Draw();
        // hist_vz
        canvas2->cd(4);
        hist_vz2->GetXaxis()->SetTitle("vz (mm)");
        hist_vz2->GetXaxis()->SetTitleSize(0.05);
        hist_vz2->GetYaxis()->SetTitle("#events");
        hist_vz2->GetYaxis()->SetTitleSize(0.05);
        hist_vz2->Draw();
        // SAVE
        canvas2->Print("../output/rec_particle.pdf");
        delete hist_p2; delete hist_theta2; delete hist_phi2; delete hist_vz2;
        delete canvas2;
    } 
    else {
        std::cout << " ***** please provide a filename..." << std::endl;
    }
    return 0;
}


