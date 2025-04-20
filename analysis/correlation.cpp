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

int layer2number(int digit) {
	if      (digit == 11) {
		return 1;
	} 
	else if (digit == 21) {
		return 2;
	} 
	else if (digit == 22) {
		return 3;
	} 
	else if (digit == 31) {
		return 4;
	} 
	else if (digit == 32) {
		return 5;
	} 
	else if (digit == 41) {
		return 6;
	} 
	else if (digit == 42) {
		return 7;
	} 
	else if (digit == 51) {
		return 8;
	} else {
		return -1; // not a layer
	}
}

int main(int argc, char const *argv[]){

	if (argc < 1) {
		printf("Please, provide a filename...");
		return 1;
	}
	
	std::vector<std::string> Data     = {"leadingEdgeTime", "timeOverThreshold", "ADC", "ped"};
	std::vector<int>         NbBin    = {20, 20, 100, 100};
	std::vector<double>      Start    = {0, 0, 0, 0};
	std::vector<double>      End      = {1000, 1000, 4000, 1000};
	std::vector<std::string> Title    = {"time", "tot", "adc", "ped"};
	double samplingTime = 1.0;
	
	const char * filename = argv[1];
	hipo::reader  reader(filename);
	hipo::banklist banklist = reader.getBanks({"AHDC::adc"});
	long unsigned int nEvent = 0;
	
	int n = Data.size();
	std::vector<std::vector<TH2D*>> corr(n, std::vector<TH2D*>(n, nullptr));
	std::vector<TH1D*> hist1d(n, nullptr);
	TH2D* occ = new TH2D("occupancy 11 (%)", "occupancy11", 99, 1.0, 100.0, 8, 1.0, 9.0);
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
		hist1d[i] = new TH1D(Title[i].c_str(), Title[i].c_str(), NbBin[i], Start[i], End[i]);
		hist1d[i]->GetXaxis()->SetTitle(Title[i].c_str());
		hist1d[i]->GetXaxis()->SetTitleSize(0.05);
		hist1d[i]->GetYaxis()->SetTitle("count");
		hist1d[i]->GetYaxis()->SetTitleSize(0.05);
	}
	
	// loop over events
	while( reader.next(banklist)){
		//printf(" ======= EVENT %ld =========\n", nEvent);
		//if (nEvent >= 20000) { break;} // process only 20k events
		if (nEvent % 10000 == 0) { printf("Start event %ld ... \n", nEvent);}
		for(int col = 0; col < banklist[0].getRows(); col++){ // loop over columns of the bankname
			int layer   = banklist[0].getInt("layer", col);
			int wire    = banklist[0].getInt("component", col);
			double time = banklist[0].getFloat("leadingEdgeTime", col)/samplingTime;
                        double tot  = banklist[0].getFloat("timeOverThreshold", col)/samplingTime;
                        double adc  = banklist[0].getInt("ADC", col); // expected adcMax without adcOffset
                        double ped  = banklist[0].getInt("ped", col); // adcOffset
			std::vector<double> value = {time, tot, adc, ped}; 
			
			if ((time >= 200) && (time <= 500) && (tot >= 350) && (tot <= 600) && (adc >= 0) && (adc <= 4095)&& (ped >= 180) && (ped <= 360) ) {	
				occ->Fill(wire, layer2number(layer));
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < i; j++) {
						corr[i][j]->Fill(value[i], value[j]);
					}
					hist1d[i]->Fill(value[i]);
				}
			}
		}
		nEvent++;
	}
	TH2D occ2 = (100.0/nEvent)*(*occ);
	occ2.SetTitle("occupancy (%)");
	occ2.GetXaxis()->SetTitle("wire number");
	occ2.GetXaxis()->SetTitleSize(0.05);
	occ2.GetYaxis()->SetTitle("layer number");
	occ2.GetYaxis()->SetTitleSize(0.05);
	TCanvas* canvas1 = new TCanvas("all","all hist 2d",1800, 990);
	//gStyle->SetOptStat("nemruo");
	gStyle->SetOptStat("");
	canvas1->Divide(3,3);
	canvas1->Update();
	TFile *f = new TFile("../output/correlation.root", "RECREATE");
	int num = 0;	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			num++;
			canvas1->cd(num);
			corr[i][j]->Draw("COLZ");
			canvas1->Update();
			corr[i][j]->Write(TString::Format("%s_%s", Title[i].c_str(), Title[j].c_str()).Data());
		}
		hist1d[i]->Write();
	}
	canvas1->cd(7);
	occ2.Draw("COLZ");
	canvas1->Update();
	canvas1->Print("../output/correlation.pdf");
	canvas1->Write();
	occ->Write("occupancy");
	occ2.Write("occupancy normalized");
	f->Close();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			delete corr[i][j];
		}
		delete hist1d[i];
	}
	delete occ;
	
	printf("nEvent : %ld\n", nEvent);
	delete canvas1;
}
