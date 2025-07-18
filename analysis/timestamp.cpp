/***********************************************
 * Sudy the effect of the timestamp on the time	
 *
 * @author Felix Touchte Codjo
 * @date July 18, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <bitset>

#include <vector>
#include <string>

#include "reader.h"
#include "futils.h"
#include "timeOffsets.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"

double fineTimestampCorrection(long timestamp, double dream_clock, bool verbose = false) {
	std::string binaryTimestamp = std::bitset<64>(timestamp).to_string();
	int fineTimestamp = (int) std::bitset<64>(binaryTimestamp.substr(64-3,3)).to_ulong();
	if (verbose) {
		printf("timestamp       : %ld\n", timestamp);
		printf("binaryTimestamp : %s\n", binaryTimestamp.c_str());
		printf("fineTimestamp   : %d\n", fineTimestamp);
		printf("timeCorrection  : %lf\n", (fineTimestamp + 0.5)*dream_clock);
	}
	return (fineTimestamp + 0.5)*dream_clock;
}

int main(int argc, char const *argv[]) {
	// test bitset library
	/*std::bitset<8> b1(10);
	std::bitset<8> b2(std::string("10000"));
	std::bitset<8> b3(0b11);
	printf("b1 : %s \n", b1.to_string().c_str());
	printf("b2 : %ld \n", b2.to_ulong());
	printf("b3 : %s \n", b3.to_string().c_str());
	printf("b2.count : %d \n", b2.size());
	fineTimestampCorrection(275540503, 8);
	return 0;*/
	if (argc > 1){
		hipo::reader  reader(argv[1]);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
		hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
		hipo::bank  eventBank(factory.getSchema("REC::Event"));
		hipo::event event;
		long unsigned int nevents =0;
		// Histograms
		TH1D* H1_leadingEdgeTime = new TH1D("leadingEdgeTime", "leadingEdgeTime (ns)", 100, 0, 1000); 
		TH1D* H1_currentTime     = new TH1D("currentTime", "leadingEdgeTime - t0 (ns)", 100, 0, 1000); 
		TH1D* H1_t0Calibration	 = new TH1D("t0_calibration", "t0_calibration (ns)", 100, 0, 400); 
		TH1D* H1_timestampCorrection  = new TH1D("timestampCorrection", "timestampCorrection (ns)", 100, 0, 60); 
		TH1D* H1_triggerTime	 = new TH1D("triggerTime", "triggerTime (ns)", 100, -1000, 400); 
		TH1D* H1_timeCorrected   = new TH1D("timeCorrected", "timeCorrected (ns)", 100, 0, 1000); 
		
		// Loop over events
		while( reader.next()){
			nevents++;
			if (nevents % 10000 == 0) { 
				printf("Start event %ld\n", nevents);
			}
			reader.read(event);
			event.getStructure(adcBank);
			event.getStructure(wfBank);
			event.getStructure(eventBank);
			for (int i = 0; i < adcBank.getRows(); i++) {
				int sector    = adcBank.getInt("sector", i);
				int layer     = adcBank.getInt("layer", i);
				int component = adcBank.getInt("component", i);
				double t0 = -9999; // fancy value
				for (int row = 0; row < timeOffsets::entries; row++) {
					if ((sector == 1) && (layer == timeOffsets::layer.at(row)) && (component == timeOffsets::component.at(row))) {
						t0 = timeOffsets::t0.at(row);
						break;
					}
				}

				double leadingEdgeTime  = adcBank.getFloat("leadingEdgeTime",i);
				double triggerTime 	= eventBank.getFloat("startTime",0);
				long timestamp 		= wfBank.getLong("timestamp",i);
				double timestampCorrection = fineTimestampCorrection(timestamp, 8, true);
				double timeCorrected       = leadingEdgeTime - t0 - triggerTime - timestampCorrection;
				H1_leadingEdgeTime->Fill(leadingEdgeTime);
				H1_t0Calibration->Fill(t0); 
				H1_currentTime->Fill(leadingEdgeTime-t0); 
				H1_triggerTime->Fill(triggerTime);
				H1_timestampCorrection->Fill(timestampCorrection);
				H1_timeCorrected->Fill(timeCorrected);
			}
		}
		printf("nevents    : %ld \n", nevents);
		// output
		const char * output = "./output/timestamp_study.root";
		TFile *f = new TFile(output, "RECREATE");
		H1_leadingEdgeTime->Write("leadingEdgeTime");
		H1_t0Calibration->Write("t0_calibration"); 
		H1_currentTime->Write("currentTime"); 
		H1_timestampCorrection->Write("timestampCorrection");
		H1_triggerTime->Write("triggerTime");
		H1_timeCorrected->Write("timeCorrected");
		f->Close();
		printf("File created : %s\n", output);
	} 
	else {
		std::cout << " ***** please provide a filename..." << std::endl;
	}

	return 0;
}


