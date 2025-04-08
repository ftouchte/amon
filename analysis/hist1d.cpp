/****************************************************
 * Generate a 1D histogram from bankname, quantity
 * to histogram, lower and upper limit
 *
 * @author Felix Touchte Codjo
 * @date March 29, 2025
 * *************************************************/

#include "reader.h"

#include <string>
#include <cstdio>

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TString.h"


int main(int argc, char const *argv[]){
	
	if (argc >= 7) { 
		// open file and read bank
		const char* filename = argv[1];
		const char* bankname = argv[2];
		const char* attribut_name = argv[3];
		const char* type = argv[4];
		int Nbins = std::atoi(argv[5]);
		double xmin = std::atof(argv[6]);
		double xmax = std::atof(argv[7]);
		
		hipo::reader  reader(filename);
		hipo::banklist banklist = reader.getBanks({bankname});
		long unsigned int nEvent = 0;
		
		char buffer[50];
		sprintf(buffer, "hist1d_%s", attribut_name);
		TH1D* hist1d = new TH1D(buffer, buffer, Nbins, xmin, xmax);	
		// loop over events
		while( reader.next(banklist)){
			//printf(" ======= EVENT %ld =========\n", nEvent);
			if (nEvent > 20000) { break;} // process only 20k events
			for(int col = 0; col < banklist[0].getRows(); col++){ // loop over columns of the bankname
				if (std::string(type) == "-f") {
					double value = banklist[0].getFloat(attribut_name, col)/50.0;
					hist1d->Fill(value);
				}
				else if (std::string(type) == "-i") {
					double value = banklist[0].getInt(attribut_name, col);
					hist1d->Fill(value);
				}
				else {
					// do nothing
				}
			}
			nEvent++;
		}
		TCanvas* canvas1 = new TCanvas("c1","c1 title",1300, 800);
		gStyle->SetOptStat("nemruo");
		hist1d->GetXaxis()->SetTitle(attribut_name);
		hist1d->GetXaxis()->SetTitleSize(0.05);
		hist1d->GetYaxis()->SetTitle("count");
		hist1d->GetYaxis()->SetTitleSize(0.05);
		hist1d->Draw();
		char buffer2[50];
		sprintf(buffer2, "../output/%s_%s.pdf", bankname, attribut_name);
		canvas1->Print(buffer2);
		delete hist1d;
		delete canvas1;
	}
	else { 
		printf("Please, all fields are mandatory...\n");
		printf("Usage :\n");
		printf("   ./hist1d filename bankname attribut type Nbins lower_value upper_value\n");
		printf("   e.g /hist1d file.hipo AHDC::adc time -f 100 0.0 100.0\n");
		return 0;
	}
}
