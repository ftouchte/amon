/***********************************************
 * Test amon/ressources/load_ccdb	
 *
 * @author Felix Touchte Codjo
 * @date July 10, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>

#include "reader.h"
#include "futils.h"
#include "timeOffsets1.h"
#include "timeOffsets2.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"


int main(int argc, char const *argv[]) {

	// Understand timeOffsets.h or load_ccdb
	TH1D* H1_t01  = new TH1D("hist_t0_1", "t0 (1)", 100, 0, 800);
	TH1D* H1_t02  = new TH1D("hist_t0_2", "t0 (2)", 100, 0, 800);
	TGraph* gr1   = new TGraph(576);
	TGraph* gr2   = new TGraph(576);
	gr1->SetTitle("t0 vs wire (1)");
	gr1->GetXaxis()->SetTitle("wire");
	gr1->GetYaxis()->SetTitle("t0");
	gr2->SetTitle("t0 vs wire (2)");
	gr2->GetXaxis()->SetTitle("wire");
	gr2->GetYaxis()->SetTitle("t0");
	for (int row = 0; row < timeOffsets1::entries; row++) {
		H1_t01->Fill(timeOffsets1::t0.at(row));
		gr1->SetPoint(row, row, timeOffsets1::t0.at(row));
	}
	for (int row = 0; row < timeOffsets2::entries; row++) {
		H1_t02->Fill(timeOffsets2::t0.at(row));
		gr2->SetPoint(row, row, timeOffsets2::t0.at(row));
	}

	// analyis of a hipo file
	if (argc < 2){
		printf("no hipo file provided; exit...\n");
		return -1;
	}
	hipo::reader  reader(argv[1]);
	hipo::dictionary factory;
	reader.readDictionary(factory);
	hipo::bank  ahdc_adc(factory.getSchema("AHDC::adc"));
	hipo::event event;
	long unsigned int nevents =0;
	// Histograms
	// electron
	TH1D* H1_leadingEdgeTime    = new TH1D("leadingEdgeTime", "leadingEdgeTime", 100, -500, 1000); 
	TH1D* H1_time    = new TH1D("time", "time : leadingEdgeTime - t0", 100, -500, 1500); 

	// Loop over events
	while( reader.next()){
		nevents++;
		//if (nevents > 10000) break;
		if (nevents % 100000 == 0) { printf("Start event %ld\n", nevents);}
		reader.read(event);
		event.getStructure(ahdc_adc);
		// loop over AHDC::adc rows
		for (int i = 0; i < ahdc_adc.getRows(); i++) {
			int sector    = ahdc_adc.getInt("sector", i);
			int layer     = ahdc_adc.getInt("layer", i);
			int component = ahdc_adc.getInt("component", i);
			double t0 = 0;
			for (int row = 0; row < timeOffsets1::entries; row++) {
				if ((sector == 1) && (layer == timeOffsets2::layer.at(row)) && (component == timeOffsets2::component.at(row))) {
					t0 = timeOffsets2::t0.at(row);
					break;
				}
			}
			double leadingEdgeTime = ahdc_adc.getFloat("leadingEdgeTime", i);
			H1_leadingEdgeTime->Fill(leadingEdgeTime);
			H1_time->Fill(leadingEdgeTime-t0);
		}
	}
	const char * output = "./output/ahdc_timeOffsets.root";
	TFile *f = new TFile(output, "RECREATE");
	H1_t01->Write("t01");
	gr1->Write("wire_vs_t01");
	H1_t02->Write("t02");
	gr2->Write("wire_vs_t02");
	H1_leadingEdgeTime->Write("leadingEdgeTime");
	H1_time->Write("time");
	f->Close();
	
	printf("File created : %s\n", output);

	return 0;
}


