/****************************************************
 * Old : histAhdcAdc.cpp  
 * from ftouchte/alert/cpp/ahdc/analysis
 * revisited to be used on files from the AHDC decoder
 *
 * @author Felix Touchte Codjo
 * @date February 17, 2025
 * *************************************************/

#include "AhdcExtractor.h"
#include "hipo4/reader.h"

#include <string>
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TString.h"


int main(int argc, char const *argv[]){
	
	if (argc < 2) {
		std::cout << "Please, provide a filename..." << std::endl;
		return 0;
	}

	using namespace std::__cxx11; // for the to_string() method
	// open file and read bank
	const char* filename = argv[1];
	hipo::reader  r(filename);
	hipo::banklist list = r.getBanks({"AHDC::adc","AHDC::wf"});
	long unsigned int nEvent =0;
	
	// Setting parameters
	double amplitudeFractionCFA = 0.5;
	double fractionCFD = 0.3;
	int binDelayCFD = 5;

	// 1D Histograms (AHDC::adc entries from extractor output)
	double tmin = 0; // time plot limit
	double tmax = 2000;
	int bins = 200;
	//TH1D* hist1d_adcMax = new TH1D("hist1d_adcMax", "adcMax (adc)",bins,0,AhdcExtractor::ADC_LIMIT);
	TH1D* hist1d_adcMax = new TH1D("hist1d_adcMax", "adcMax (adc)",bins,0, 1000);
	TH1D* hist1d_integral = new TH1D("hist1d_integral", "integral (adc per 44 ns)",bins,0,6); // we store log10(integral) 
	TH1D* hist1d_adcOffset = new TH1D("hist1d_adcOffset", "adcOffest (adc)",bins,100,500);
	TH1D* hist1d_timeMax = new TH1D("hist1d_timeMax", "timeMax (ns)",bins,0,60);
	TH1D* hist1d_timeRiseCFA = new TH1D("hist1d_timeRiseCFA", "leadingEdgeTime",bins, 0.0, 60.0);
	TH1D* hist1d_timeOT = new TH1D("hist1d_timeOT", "timeOverThreshold",bins, 0.0, 60.0);
	TH1D* hist1d_timeCFD = new TH1D("hist1d_timeCFD", "constantFractionTimeCFD",bins,1000, 60.0);

	// loop over events
	while( r.next(list)){
		if (nEvent % 1000 == 0)
			std::cout << "***** starting event " << nEvent << std::endl;
		// loop over columns of AHDC::wf
		for(int col = 0; col < list[1].getRows(); col++){ 
			std::vector<short> samples; // waveform samples
			double my_adcMax = 0.0;
			for (int bin=0; bin < 128; bin++){
				std::string binName = "s" + to_string(bin+1);
				short thisSample = list[1].getInt(binName.c_str(),col);
				samples.push_back(thisSample);
				my_adcMax = (my_adcMax < thisSample) ? thisSample : my_adcMax;
			}
			//AhdcExtractor T(44.0,amplitudeFractionCFA,binDelayCFD,fractionCFD);
			AhdcExtractor T(1.0,amplitudeFractionCFA,binDelayCFD,fractionCFD);
			T.adcOffset = (short) (samples[0] + samples[1] + samples[2] + samples[3] + samples[4])/5;
			std::map<std::string,double> output = T.extract(samples);
			
			// Fill histograms
			double adcMax = output["adcMax"];
			double integral = output["integral"];
			double adcOffset = output["adcOffset"];
			double timeMax = output["timeMax"];
			double timeRiseCFA = output["leadingEdgeTime"];
			double timeOT = output["timeOverThreshold"];
			double timeCFD = output["constantFractionTime"];
			//hist1d_adcMax->Fill(adcMax);
			hist1d_adcMax->Fill(my_adcMax);
			hist1d_integral->Fill(integral > 1 ? log10(integral) : 0);
			hist1d_adcOffset->Fill(adcOffset);
			hist1d_timeMax->Fill(timeMax);
			hist1d_timeRiseCFA->Fill(timeRiseCFA);
			hist1d_timeCFD->Fill(timeCFD);
			hist1d_timeOT->Fill(timeOT);
		} // end loop over columns of AHDC::wf
		nEvent++;
	} // end loop over events

	// Plots
	//TCanvas* canvas1 = new TCanvas("c1","c1 title",2000,3000);
	TCanvas* canvas1 = new TCanvas("c1","c1 title",3500,2000);
	canvas1->Divide(2,4); // (nx,ny)
	gStyle->SetOptStat("nemruo");

	// adcMax
        canvas1->cd(1);
        hist1d_adcMax->GetXaxis()->SetTitle("adcMax (adc)");
        hist1d_adcMax->GetXaxis()->SetTitleSize(0.05);
        hist1d_adcMax->GetYaxis()->SetTitle("count");
        hist1d_adcMax->GetYaxis()->SetTitleSize(0.05);
        hist1d_adcMax->Draw();
	
	// integral 
        canvas1->cd(2);
	//gPad->SetLogx(1);
        hist1d_integral->GetXaxis()->SetTitle("Log10 integral (adc per 44 ns)");
        hist1d_integral->GetXaxis()->SetTitleSize(0.05);
        hist1d_integral->GetYaxis()->SetTitle("count");
        hist1d_integral->GetYaxis()->SetTitleSize(0.05);
        hist1d_integral->Draw();

	// adcOffset
        canvas1->cd(3);
        hist1d_adcOffset->GetXaxis()->SetTitle("adcOffset (adc)");
        hist1d_adcOffset->GetXaxis()->SetTitleSize(0.05);
        hist1d_adcOffset->GetYaxis()->SetTitle("count");
        hist1d_adcOffset->GetYaxis()->SetTitleSize(0.05);
        hist1d_adcOffset->Draw();
	
	
	// timeMax
        canvas1->cd(5);
        hist1d_timeMax->GetXaxis()->SetTitle("timeMax");
        hist1d_timeMax->GetXaxis()->SetTitleSize(0.05);
        hist1d_timeMax->GetYaxis()->SetTitle("count");
        hist1d_timeMax->GetYaxis()->SetTitleSize(0.05);
        hist1d_timeMax->Draw();

	// timeRiseCFA
        canvas1->cd(6);
        hist1d_timeRiseCFA->GetXaxis()->SetTitle("leadingEdgeTime");
        hist1d_timeRiseCFA->GetXaxis()->SetTitleSize(0.05);
        hist1d_timeRiseCFA->GetYaxis()->SetTitle("count");
        hist1d_timeRiseCFA->GetYaxis()->SetTitleSize(0.05);
        hist1d_timeRiseCFA->Draw();
	
	// timeCFD
        canvas1->cd(7);
        hist1d_timeCFD->GetXaxis()->SetTitle("constantFractionTime");
        hist1d_timeCFD->GetXaxis()->SetTitleSize(0.05);
        hist1d_timeCFD->GetYaxis()->SetTitle("count");
        hist1d_timeCFD->GetYaxis()->SetTitleSize(0.05);
        hist1d_timeCFD->Draw();
	
	// timeOT
        canvas1->cd(8);
        hist1d_timeOT->GetXaxis()->SetTitle("timeOverThreshold");
        hist1d_timeOT->GetXaxis()->SetTitleSize(0.05);
        hist1d_timeOT->GetYaxis()->SetTitle("count");
        hist1d_timeOT->GetYaxis()->SetTitleSize(0.05);
        hist1d_timeOT->Draw();
	
	// output
	canvas1->Print(TString::Format("histAhdcAdc.pdf"));
	delete canvas1;
	/*********************
	 * canvas 2
	 * ******************/
	TCanvas* canvas2 = new TCanvas("c2","c2 title",3500,2000);
        canvas2->Divide(2,2); // (nx,ny)
        gStyle->SetOptStat("nemruo");
	canvas2->cd(1);
	hist1d_timeMax->Draw();
	canvas2->cd(2);
	hist1d_timeRiseCFA->Draw();
	canvas2->cd(3);
	hist1d_timeCFD->Draw();
	canvas2->cd(4);
	hist1d_timeOT->Draw();
        canvas2->Print(TString::Format("histAhdcAdc2.pdf"));
	delete canvas2;

	// Delete
	delete hist1d_adcMax; 
	delete hist1d_integral;
	delete hist1d_adcOffset;
	delete hist1d_timeMax;
	delete hist1d_timeRiseCFA;
	delete hist1d_timeOT;
	delete hist1d_timeCFD;
}
