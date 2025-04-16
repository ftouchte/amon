/****************************************************
 * Study of the correlation between the decoding
 * parameters
 *
 * @author Felix Touchte Codjo
 * @date April 15, 2025
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

	if (argc < 1) {
		printf("Please, provide a filename...");
		return 1;
	}
	
	std::vector<std::string> Data     = {"ADC", "ped", "leadingEdgeTime", "timeMax", "timeOverThreshold"};
	std::vector<int>         NbBin    = {100, 100, 20, 20, 20};
	std::vector<double>      Start    = {0.0, 0.0, 0.0, 0.0, 0.0};
	std::vector<double>      End      = {4095.0, 1000.0, 20.0, 20.0, 20.0};
	std::vector<std::string> Title    = {"adcMax", "adcOffset", "leadingEdgeTime", "timeMax", "timeOverThreshold"};
	double samplingTime = 50.0;
	
	const char * filename = argv[1];
	hipo::reader  reader(filename);
	hipo::banklist banklist = reader.getBanks({"AHDC::adc"});
	long unsigned int nEvent = 0;
	
	int n = Data.size();
	std::vector<std::vector<TH2D*>> corr(n, std::vector<TH2D*>(n, nullptr));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			char xtitle[80];
			char ytitle[80];
			sprintf(xtitle, "%s %d %d", Title[i].c_str(), i, j);
			sprintf(ytitle, "%s %d %d", Title[j].c_str(), i, j);
			corr[i][j] = new TH2D(xtitle, ytitle, NbBin[i], Start[i], End[i], NbBin[j], Start[j], End[j]);
			// title
			corr[i][j]->SetTitle(TString::Format("%s vs %s", xtitle, ytitle).Data());
			corr[i][j]->GetXaxis()->SetTitle(xtitle);
			corr[i][j]->GetXaxis()->SetTitleSize(0.05);
			corr[i][j]->GetYaxis()->SetTitle(ytitle);
			corr[i][j]->GetYaxis()->SetTitleSize(0.05);
		}
	}
	
	// loop over events
	while( reader.next(banklist)){
		//printf(" ======= EVENT %ld =========\n", nEvent);
		if (nEvent >= 20000) { break;} // process only 20k events
		for(int col = 0; col < banklist[0].getRows(); col++){ // loop over columns of the bankname
			double timeMax = banklist[0].getFloat("time", col)/samplingTime;
                        double leadingEdgeTime = banklist[0].getFloat("leadingEdgeTime", col)/samplingTime;
                        double timeOverThreshold = banklist[0].getFloat("timeOverThreshold", col)/samplingTime;
                        double constantFractionTime = banklist[0].getFloat("constantFractionTime", col)/samplingTime;
                        double adcMax = banklist[0].getInt("ADC", col); // expected adcMax without adcOffset
                        double adcOffset = banklist[0].getInt("ped", col);
			std::vector<double> value = {adcMax, adcOffset, leadingEdgeTime, timeMax, timeOverThreshold, constantFractionTime}; 
			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < i; j++) {
					corr[i][j]->Fill(value[i], value[j]);
				}
			}
		}
		nEvent++;
	}
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1800, 990);
	//gStyle->SetOptStat("nemruo");
	gStyle->SetOptStat("");
	int n_hits = n*(n-1)/2;
	//canvas1->Divide(3, n_hits/3);
	canvas1->Divide(3,3);
	canvas1->Update();
	int num = 0;	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			num++;
			if (num > 9) break;
			printf("num : %d\n", num);
			canvas1->cd(num);
			corr[i][j]->Draw("COLZ");
			canvas1->Update();
		}
	}
	canvas1->Print("../output/correlation.pdf");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			delete corr[i][j];
		}
	}
	printf("nEvent : %ld\n", nEvent);
	delete canvas1;
}
