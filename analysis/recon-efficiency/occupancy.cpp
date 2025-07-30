/***********************************************
 * Occupancy study (raw hit cuts vs wfType)
 *
 * A given file (see filename var.)
 * - Plot the occupancy in percentage (with respect to the number of event)
 *    * Without cuts
 *    * with my raw hit cuts on (time, tot, adc, ped... see ccdb)
 *    * with Noemie wfType
 * - Plot the disrition of time, tot, adc and ped according to the above cuts
 *
 * Remarks: test on run 22435/00003
 * - wfType is a good compromise (6 %); rawCuts (< 4%)
 * - wfType <= 1 keep more track compared the rawCuts at low adc
 * - wfType <= 1 keep less track compared the rawCuts at large adc
 * - more time cut can be added to complete wfType <= 1
 * - starting wfType 2, we have a default value of ped that is in the range of my rawCuts; maybe we keep bad events
 *
 * @author Felix Touchte Codjo
 * @date July 29, 2025
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
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"

//Waveform types:
//0 is good, 
//-1 is invalid,
//1 is saturating,
//2 has too short of a baseline,
//3 is late and only has a rising edge,
//4 is a trailing edge from a previous wf,
//5 has low ADC ("flat")

int layer2number(int digit);
void progressBar(int state);

int main(int argc, char const *argv[]) {
        const char * filename = "/home/touchte-codjo/Desktop/hipofiles/occupancy/new_clas_022435.00003.hipo";
		hipo::reader  reader(filename);
		hipo::dictionary factory;
		reader.readDictionary(factory);
		hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
		hipo::event event;
		long unsigned int nevents =0;
		// Histograms
        std::vector<TH1D*> wfType_histos;
        wfType_histos.push_back(new TH1D("occ_wfType0", "occupancy (wfType == 0)", 576, 1, 577));
        wfType_histos.push_back(new TH1D("occ_wfType1", "occupancy (wfType <= 1)", 576, 1, 577));
        wfType_histos.push_back(new TH1D("occ_wfType2", "occupancy (wfType <= 2)", 576, 1, 577));
        wfType_histos.push_back(new TH1D("occ_wfType3", "occupancy (wfType <= 3)", 576, 1, 577));
        wfType_histos.push_back(new TH1D("occ_wfType4", "occupancy (wfType <= 4)", 576, 1, 577));
        wfType_histos.push_back(new TH1D("occ_wfType5", "occupancy (wfType <= 5)", 576, 1, 577));
		TH1D* rawCuts_histo  = new TH1D("occ_rawCuts", "occupancy (rawCuts)", 576, 1, 577); 
		TH1D* noCuts_histo   = new TH1D("occ_noCuts", "occupancy (noCuts)", 576, 1, 577);
        // Decoding quantities
        std::vector<TH1D*> TIME;
        TIME.push_back(new TH1D("time_wftype0", "time (wfType == 0)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_wftype1", "time (wfType <= 1)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_wftype2", "time (wfType <= 2)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_wftype3", "time (wfType <= 3)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_wftype4", "time (wfType <= 4)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_wftype5", "time (wfType <= 5)", 100, 0, 1000));
        TIME.push_back(new TH1D("time_rawCuts", "time (rawCuts)", 100, 0, 1000)); // raw cuts
        TIME.push_back(new TH1D("time_noCuts" , "time (noCuts)" , 100, 0, 1000)); // no cuts
        std::vector<TH1D*> TOT;
        TOT.push_back(new TH1D("tot_wftype0", "tot (wfType == 0)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_wftype1", "tot (wfType <= 1)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_wftype2", "tot (wfType <= 2)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_wftype3", "tot (wfType <= 3)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_wftype4", "tot (wfType <= 4)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_wftype5", "tot (wfType <= 5)", 100, 0, 1000));
        TOT.push_back(new TH1D("tot_rawCuts", "tot (rawCuts)", 100, 0, 1000)); // raw cuts
        TOT.push_back(new TH1D("tot_noCuts" , "tot (noCuts)" , 100, 0, 1000)); // no cuts
        std::vector<TH1D*> ADC;
        ADC.push_back(new TH1D("adc_wftype0", "adc (wfType == 0)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_wftype1", "adc (wfType <= 1)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_wftype2", "adc (wfType <= 2)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_wftype3", "adc (wfType <= 3)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_wftype4", "adc (wfType <= 4)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_wftype5", "adc (wfType <= 5)", 100, 0, 3700));
        ADC.push_back(new TH1D("adc_rawCuts", "adc (rawCuts)", 100, 0, 3700)); // raw cuts
        ADC.push_back(new TH1D("adc_noCuts" , "adc (noCuts)" , 100, 0, 3700)); // no cuts
        std::vector<TH1D*> PED;
        PED.push_back(new TH1D("ped_wftype0", "ped (wfType == 0)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_wftype1", "ped (wfType <= 1)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_wftype2", "ped (wfType <= 2)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_wftype3", "ped (wfType <= 3)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_wftype4", "ped (wfType <= 4)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_wftype5", "ped (wfType <= 5)", 100, 0, 1000));
        PED.push_back(new TH1D("ped_rawCuts", "ped (rawCuts)", 100, 0, 1000)); // raw cuts
        PED.push_back(new TH1D("ped_noCuts" , "ped (noCuts)" , 100, 0, 1000)); // no cuts

		// Loop over events
		while( reader.next()){
			nevents++;
            if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
                progressBar(100.0*nevents/reader.getEntries());
            }
			reader.read(event);
			event.getStructure(adcBank);
			for (int i = 0; i < adcBank.getRows(); i++) {
			    // sufficient for wfType study
                int wfType = adcBank.getInt("wfType", i);
                int layer  = adcBank.getInt("layer", i);
                int component = adcBank.getInt("component", i);
                // relevant raw hit cuts
                double time = adcBank.getFloat("leadingEdgeTime", i);
                double tot  = adcBank.getFloat("timeOverThreshold", i);
                double adc  = adcBank.getInt("ADC", i); // expected adcMax without adcOffset
                double ped  = adcBank.getFloat("ped", i);
                // Fill histograms
                // layer 1 --> 47 wires 
                // layer 2 --> 56 wires
                // layer 3 --> 56 wires 
                // layer 4 --> 72 wires
                // layer 5 --> 72 wires
                // layer 6 --> 87 wires
                // layer 7 --> 87 wires
                // layer 8 --> 99 wires
                int layerNumber = layer2number(layer);
                int wire = -1;
                if (layerNumber == 1) {
                    wire = component;
                }
                else if (layerNumber == 2) {
                    wire = 47 + component;
                }
                else if (layerNumber == 3) {
                    wire = 47 + 56 + component;
                }
                else if (layerNumber == 4) {
                    wire = 47 + 56 + 56 + component;
                }
                else if (layerNumber == 5) {
                    wire = 47 + 56 + 56 + 72 + component;
                }
                else if (layerNumber == 6) {
                    wire = 47 + 56 + 56 + 72 + 72 + component;
                }
                else if (layerNumber == 7) {
                    wire = 47 + 56 + 56 + 72 + 72 + 87 + component;
                }
                else if (layerNumber == 8) {
                    wire = 47 + 56 + 56 + 72 + 72 + 87 + 87 + component;
                }
                else {
                    wire = -11; // it is getting worst actually
                }
                // wfType cuts
                if (wfType == 0) {
                    wfType_histos[0]->Fill(wire);
                    TIME[0]->Fill(time);
                    TOT[0]->Fill(tot);
                    ADC[0]->Fill(adc);
                    PED[0]->Fill(ped);
                }
                for (int j : std::vector<int>({1, 2, 3, 4, 5})){
                    if (wfType <= j) {
                        wfType_histos[j]->Fill(wire);
                        TIME[j]->Fill(time);
                        TOT[j]->Fill(tot);
                        ADC[j]->Fill(adc);
                        PED[j]->Fill(ped);
                    }
                }
                // raw hit cut histo
                if ((time >= 0) && (time <= 300) && (tot >= 300) && (tot <= 600) && (ped >= 180) && (ped <= 360) && (adc >= 0) && (adc <= 4095)) {
                    rawCuts_histo->Fill(wire);
                    TIME[6]->Fill(time);
                    TOT[6]->Fill(tot);
                    ADC[6]->Fill(adc);
                    PED[6]->Fill(ped);
                }
                // no cut histo
                noCuts_histo->Fill(wire);
                TIME[7]->Fill(time);
                TOT[7]->Fill(tot);
                ADC[7]->Fill(adc);
                PED[7]->Fill(ped);
            }
		}
		printf("nevents    : %ld \n", nevents);
        // Renormalize over the number of events
        for (auto ptr_histo : wfType_histos) {
            ptr_histo->Scale(100.0/nevents);
        }
        rawCuts_histo->Scale(100.0/nevents);
        noCuts_histo->Scale(100.0/nevents);
        // Label axis
        for (auto ptr : wfType_histos) {
            ptr->GetXaxis()->SetTitle("wire");
            ptr->GetXaxis()->SetTitleSize(0.05);;
            ptr->GetYaxis()->SetTitle("count [%]");
            ptr->GetYaxis()->SetTitleSize(0.05);;
        }
        rawCuts_histo->GetXaxis()->SetTitle("wire");
        rawCuts_histo->GetXaxis()->SetTitleSize(0.05);;
        rawCuts_histo->GetYaxis()->SetTitle("count [%]");
        rawCuts_histo->GetYaxis()->SetTitleSize(0.05);;
        noCuts_histo->GetXaxis()->SetTitle("wire");
        noCuts_histo->GetXaxis()->SetTitleSize(0.05);;
        noCuts_histo->GetYaxis()->SetTitle("count [%]");
        noCuts_histo->GetYaxis()->SetTitleSize(0.05);
        for (auto ptr : TIME) {
            ptr->GetXaxis()->SetTitle("time (ns)");
            ptr->GetXaxis()->SetTitleSize(0.05);;
            ptr->GetYaxis()->SetTitle("count");
            ptr->GetYaxis()->SetTitleSize(0.05);;
        }
        for (auto ptr : TOT) {
            ptr->GetXaxis()->SetTitle("tot (ns)");
            ptr->GetXaxis()->SetTitleSize(0.05);;
            ptr->GetYaxis()->SetTitle("count");
            ptr->GetYaxis()->SetTitleSize(0.05);;
        }
        for (auto ptr : ADC) {
            ptr->GetXaxis()->SetTitle("adc");
            ptr->GetXaxis()->SetTitleSize(0.05);;
            ptr->GetYaxis()->SetTitle("count");
            ptr->GetYaxis()->SetTitleSize(0.05);;
        }
        for (auto ptr : PED) {
            ptr->GetXaxis()->SetTitle("ped");
            ptr->GetXaxis()->SetTitleSize(0.05);;
            ptr->GetYaxis()->SetTitle("count");
            ptr->GetYaxis()->SetTitleSize(0.05);;
        }
		// output
		const char * output = "./occupancy.root";
		TFile *f = new TFile(output, "RECREATE");
        TDirectory *occupancy = f->mkdir("occupancy");
        TDirectory *time = f->mkdir("time");
        TDirectory *tot = f->mkdir("tot");
        TDirectory *adc = f->mkdir("adc");
        TDirectory *ped = f->mkdir("ped");
        occupancy->cd();
		wfType_histos[0]->Write("wfType==0");
		wfType_histos[1]->Write("wfType<=1");
		wfType_histos[2]->Write("wfType<=2");
		wfType_histos[3]->Write("wfType<=3");
		wfType_histos[4]->Write("wfType<=4");
		wfType_histos[5]->Write("wfType<=5");
        rawCuts_histo->Write("rawCuts");
        noCuts_histo->Write("noCuts");
        time->cd();
		TIME[0]->Write("wfType==0");
		TIME[1]->Write("wfType<=1");
		TIME[2]->Write("wfType<=2");
		TIME[3]->Write("wfType<=3");
		TIME[4]->Write("wfType<=4");
		TIME[5]->Write("wfType<=5");
		TIME[6]->Write("rawCuts");
		TIME[7]->Write("noCuts");
        tot->cd();
		TOT[0]->Write("wfType==0");
		TOT[1]->Write("wfType<=1");
		TOT[2]->Write("wfType<=2");
		TOT[3]->Write("wfType<=3");
		TOT[4]->Write("wfType<=4");
		TOT[5]->Write("wfType<=5");
		TOT[6]->Write("rawCuts");
		TOT[7]->Write("noCuts");
        adc->cd();
		ADC[0]->Write("wfType==0");
		ADC[1]->Write("wfType<=1");
		ADC[2]->Write("wfType<=2");
		ADC[3]->Write("wfType<=3");
		ADC[4]->Write("wfType<=4");
		ADC[5]->Write("wfType<=5");
		ADC[6]->Write("rawCuts");
		ADC[7]->Write("noCuts");
        ped->cd();
		PED[0]->Write("wfType==0");
		PED[1]->Write("wfType<=1");
		PED[2]->Write("wfType<=2");
		PED[3]->Write("wfType<=3");
		PED[4]->Write("wfType<=4");
		PED[5]->Write("wfType<=5");
		PED[6]->Write("rawCuts");
		PED[7]->Write("noCuts");
		f->Close();
		printf("File created : %s\n", output);
        // Memory deallocation
        for (auto ptr : wfType_histos) {
            delete ptr;
        }
        delete rawCuts_histo;
        delete noCuts_histo;
        for (auto ptr : TIME) {
            delete ptr;
        }
        for (auto ptr : TOT) {
            delete ptr;
        }
        for (auto ptr : ADC) {
            delete ptr;
        }
        for (auto ptr : PED) {
            delete ptr;
        }

	return 0;
}

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

void progressBar(int state) { // state is a number between 0 and 100
    if (state > 100) { return;}
    printf("\rProgress \033[32m\[");
    for (int i = 0; i <= state; i++) {
        printf("#");
    }
    printf("\033[0m");
    for (int i = state+1; i < 100; i++) {
        printf(".");
    }
    if (state == 100) {
        printf("\033[32m] \033[1m %d %%\033[0m\n", state);
    } else {
        printf("] %d %%", state);
    }
    fflush(stdout);
}
