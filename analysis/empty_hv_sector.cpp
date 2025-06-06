/***********************************************
 * Study AHDC HV sector occupancy
 *
 * @author Felix Touchte Codjo
 * @date May 26, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>

#include "reader.h"
#include "futils.h"
#include "AhdcMapping.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TString.h"

namespace futils {
	void Get_HV_sector(int sector, int layer, int component, int & crate, int & slot, int & channel, int & hv, int & sub_hv);
	int layer2number(int digit);
}

int main(int argc, char const *argv[]) {
	if (argc > 1){
		printf("File name : %s", argv[1]);
		hipo::reader  reader(argv[1]);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
		hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
		hipo::event event;
		long unsigned int nevents = 0;
		// Histograms
		TH2D* H2_occupancy = new TH2D("occupancy", "occupancy", 99, 1, 100, 8, 1, 9); 
		TH2D* H2_adcMax    = new TH2D("adcMax", "adcMax", 99, 1, 100, 8, 1, 9); 
		TH2D* H2_time      = new TH2D("time", "time", 99, 1, 100, 8, 1, 9); 
		TH2D* H2_tot       = new TH2D("tot", "tot", 99, 1, 100, 8, 1, 9); 
		TH1D* H1_avg_wf    = new TH1D("wf_all", "wf_all", 20, 0, 20);
		TH1D* H1_avg_wf0    = new TH1D("wf_0", "wf_0", 20, 0, 20); // for a given wire
		std::vector<TH1D*> VectorH1_avg_wf(27, nullptr); // one per hv sector
		for (int hv = 0; hv < 9; hv++) {
			for (int shv = 0; shv < 3; shv++) {
				int index = 3*hv + shv;
				VectorH1_avg_wf[index] = new TH1D(TString::Format("HV_%d-%d", hv+1, shv+1).Data(),TString::Format("HV_%d-%d", hv+1, shv+1).Data(), 20, 0, 20);
			}
		}

		// Loop over events
		while( reader.next()){
			nevents++;
			//if (nevents > 10000) break;
			if (nevents % 1000 == 0) { printf("Start event %ld\n", nevents);}
			reader.read(event);
			event.getStructure(adcBank);
			event.getStructure(wfBank);
			for (int i = 0; i < adcBank.getRows(); i++) {
				int layer  = adcBank.getInt("layer", i);
				int comp   = adcBank.getInt("component", i);
				double adc  = adcBank.getInt("ADC", i);	
				double time = adcBank.getFloat("leadingEdgeTime", i);	
				double tot  = adcBank.getFloat("timeOverThreshold", i);
				// Fill histograms
				int layer_number = futils::layer2number(layer);
				H2_occupancy->Fill(comp, layer_number);
				H2_adcMax->Fill(comp, layer_number, adc);
				H2_time->Fill(comp, layer_number, time);
				H2_tot->Fill(comp, layer_number, tot);
				int crate, slot, channel, hv, sub_hv; 
				futils::Get_HV_sector(1, layer, comp, crate, slot, channel, hv, sub_hv);
				// Loop over samples
				for (int s = 0; s < 20; s++) {
					char sname[50];
					sprintf(sname, "s%d", s+1);
					double value = wfBank.getInt(sname, i);
					int index = 3*(hv-1) + sub_hv - 1; // here number start at 1
					VectorH1_avg_wf[index]->Fill(s, value); // per hv-sub_hv sector
					H1_avg_wf->Fill(s, value); // global
					if ((layer == 11) && (comp == 4)) {
						H1_avg_wf0->Fill(s, value);
					}
				}
			}
		}
		printf("nevents    : %ld \n", nevents);
		// output
		const char * output = "../output/hv_sector_study.root";
		TFile *f = new TFile(output, "RECREATE");
		H2_occupancy->Write("H2_occupancy");
		H2_adcMax->Write("H2_adcMax");
		H2_time->Write("H2_time");
		H2_tot->Write("H2_tot");
		(*H2_adcMax/(*H2_occupancy)).Write("H2_adcMax_avg");
		(*H2_time/(*H2_occupancy)).Write("H2_time_avg");
		(*H2_tot/(*H2_occupancy)).Write("H2_tot_avg");
		H1_avg_wf->Write("H1_avg_wf");
		H1_avg_wf0->Write("H1_avg_wf0");
		for (int hv = 0; hv < 9; hv++) {
			for (int shv = 0; shv < 3; shv++) {
				int index = 3*hv + shv;
				VectorH1_avg_wf[index]->Write(TString::Format("H1_wf_hv_%d-%d", hv+1, shv+1).Data());
			}
		}
		f->Close();
		printf("File created : %s\n", output);
	} 
	else {
		std::cout << " ***** please provide a filename..." << std::endl;
	}

	return 0;
}

int futils::layer2number(int digit) {
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

void futils::Get_HV_sector(int sector, int layer, int component, int & crate, int & slot, int & channel, int & hv, int & sub_hv) {
        AhdcMapping::GetDreamChannel(1, layer, component, crate, slot, channel);
        hv = -1;
        sub_hv = -11; 
        // read hv_notes.txt
        if (slot == 1) { 
                hv = (channel-1)/64 + 1; 
                if (hv % 2 == 0) hv--;
                else hv++;
                int num = (channel-1)%64 + 1; 
                if ((num >= 1) && (num <= 23)) {
                        sub_hv = 3; 
                }    
                else if ((num >= 24) && (num <= 45)) {
                        sub_hv = 2; 
                }    
                else if ((num >= 46) && (num <= 64)) {
                        sub_hv = 1; 
                }
        }    
        else if (slot == 4) { // hv only 9
                // read hv_notes.txt
                hv = (channel-1)/64 + 1; 
                if (hv == 2) { hv = 9;}
                if (hv != 9) {
                        if (hv % 2 == 0) hv--;
                        else hv++;
                }    
                int num = (channel-1)%64 + 1; 
                if ((num >= 1) && (num <= 23)) {
                        sub_hv = 3; 
                }
                else if ((num >= 24) && (num <= 45)) {
                        sub_hv = 2; 
                }
                else if ((num >= 46) && (num <= 64)) {
                        sub_hv = 1; 
                }
        }
}
