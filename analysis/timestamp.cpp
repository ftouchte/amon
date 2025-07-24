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
		TH1D* H1_currentTime     = new TH1D("currentTime", "leadingEdgeTime - t0 (ns)", 100, 0, 400); 
		TH1D* H1_t0Calibration	 = new TH1D("t0_calibration", "t0_calibration (ns)", 100, 150, 400); 
		TH1D* H1_timestampCorrection  = new TH1D("timestampCorrection", "FineTimestampCorrection (ns)", 100, 0, 40); 
		TH1D* H1_triggerTime	 = new TH1D("triggerTime", "triggerTime (ns)", 100, 0, 200); 
		TH1D* H1_timeCorrected   = new TH1D("timeCorrected", "timeCorrected (ns)", 100, -180, 400); 
		
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
			double triggerTime 	= eventBank.getFloat("startTime",0);
			// https://github.com/JeffersonLab/coatjava/tree/5912ac9f93b1fe5cdb7a72fbb2d370fe7c14ff48/reconstruction/eb#event-start-time
			if (abs(triggerTime + 1000) <= 1e-15*std::max(std::fabs(triggerTime), 1000.0)) {
				//printf("bad event, triggerTime : %lf\n", triggerTime);
				continue;
			}
			for (int i = 0; i < adcBank.getRows(); i++) {
				//printf("good event, triggerTime : %lf\n", triggerTime);
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
				long   timestamp 	= wfBank.getLong("timestamp",i);
				double timestampCorrection = fineTimestampCorrection(timestamp, 8);
				double timeCorrected       = leadingEdgeTime - t0 - triggerTime + timestampCorrection;
				// Raw hit cuts
				double time = leadingEdgeTime - t0;
				double tot  = adcBank.getFloat("timeOverThreshold",i);
				double ped  = adcBank.getInt("ped",i);
				double adc  = adcBank.getInt("ADC",i);
				if ((time >= 0) && (time <= 300) && (tot >= 300) && (tot <= 600) && (ped >= 180) && (ped <= 360) && (adc >= 0) && (adc <= 4095)) {
					H1_leadingEdgeTime->Fill(leadingEdgeTime);
					H1_t0Calibration->Fill(t0); 
					H1_currentTime->Fill(time); 
					H1_triggerTime->Fill(triggerTime);
					H1_timestampCorrection->Fill(timestampCorrection);
					H1_timeCorrected->Fill(timeCorrected);
				}
			}
		}
		printf("nevents    : %ld \n", nevents);
		// output
		const char * output = "./output/timestamp_study.root";
		TFile *f = new TFile(output, "RECREATE");
		H1_leadingEdgeTime->Write("leadingEdgeTime");
		H1_currentTime->Write("currentTime"); 
		H1_t0Calibration->Write("t0_calibration"); 
		H1_timestampCorrection->Write("fineTimestampCorrection");
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


