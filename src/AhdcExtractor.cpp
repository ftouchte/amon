#include "AhdcExtractor.h"


#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1.h"
#include "TGraphPolar.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TArrow.h"


std::map<std::string,double> AhdcExtractor::extract(const std::vector<short> samples){
	samplesCorr = samples;
	this->waveformCorrection();
	this->fitAverage();
	this->fitParabolic();
	this->computeTimeAtConstantFractionAmplitude();
	this->computeTimeUsingConstantFractionDiscriminator();
	//this->fineTimeStampCorrection();
	std::map<std::string,double> output;
	output["binMax"] = binMax;
	output["binOffset"] = binOffset;
	output["adcMax"] = adcMax;
	output["timeMax"] = timeMax;
	output["integral"] = integral;
	output["leadingEdgeTime"] = leadingEdgeTime;
	output["trailingEdgeTime"] = trailingEdgeTime;
	output["timeOverThreshold"] = timeOverThreshold;
	output["constantFractionTime"] = constantFractionTime;
	output["adcOffset"] = adcOffset;
	return output;

}

void AhdcExtractor::waveformCorrection(){
	binNumber = samplesCorr.size();
	binMax = 0;
	adcMax = (short) (samplesCorr[0] - adcOffset);
	integral = 0;
	for (int bin = 0; bin < binNumber; bin++){
		samplesCorr[bin] = (short) (samplesCorr[bin] - adcOffset);
		if (adcMax < samplesCorr[bin]){
			adcMax = samplesCorr[bin];
			binMax = bin;
		}
		integral += samplesCorr[bin];
	}
	/*
	 * If adcMax + adcOffset == ADC_LIMIT, that means there is saturation
	 * In that case, binMax is the middle of the first plateau
	 * This convention can be changed
	 */
	if ((short) adcMax + adcOffset == ADC_LIMIT) {
		int binMax2 = binMax;
		for (int bin = binMax; bin < binNumber; bin++){
			if (samplesCorr[bin] + adcOffset == ADC_LIMIT) {
				binMax2 = bin;
			}
			else {
				break;
			}
		}
		binMax = (binMax + binMax2)/2;
	}
	binOffset = sparseSample*binMax;
	timeMax = (binMax + binOffset)*samplingTime;
}


void AhdcExtractor::fitAverage(){
	if ((binMax - 2 >= 0) && (binMax + 2 <= binNumber - 1)){
		adcMax = 0;
		for (int bin = binMax - 2; bin <= binMax + 2; bin++){
			adcMax += samplesCorr[bin];
		}
		adcMax = adcMax/5;
	}
}

void AhdcExtractor::fitParabolic(){}

void AhdcExtractor::fineTimeStampCorrection(){}

void AhdcExtractor::computeTimeAtConstantFractionAmplitude(){
	float threshold = amplitudeFractionCFA*adcMax;
	// leadingEdgeTime
	int binRise = 0;
	for (int bin = 0; bin < binMax; bin++){
		if (samplesCorr[bin] < threshold)
			binRise = bin;  // last pass below threshold and before adcMax
	} // at this stage : binRise < leadingEdgeTime/samplingTime <= binRise + 1 // leadingEdgeTime is determined by assuming a linear fit between binRise and binRise + 1
	float slopeRise = 0;
	if (binRise + 1 <= binNumber-1)
		slopeRise = samplesCorr[binRise+1] - samplesCorr[binRise];
	float fittedBinRise = (slopeRise == 0) ? binRise : binRise + (threshold - samplesCorr[binRise])/slopeRise;
	leadingEdgeTime = (fittedBinRise + binOffset)*samplingTime; // binOffset is determined in wavefromCorrection() // must be the same for all time ? // or must be defined using fittedBinRise*sparseSample

	// trailingEdgeTime
	int binFall = binMax;
	for (int bin = binMax; bin < binNumber; bin++){
		if (samplesCorr[bin] > threshold){
				binFall = bin;
		}
		else {
				binFall = bin;
				break; // first pass below the threshold
		}
	} // at this stage : binFall - 1 <= leadingEdgeTime/samplingTime < binFall // trailingEdgeTime is determined by assuming a linear fit between binFall - 1 and binFall
	float slopeFall = 0;
	if (binFall - 1 >= 0)
		slopeFall = samplesCorr[binFall] - samplesCorr[binFall-1];
	float fittedBinFall = (slopeFall == 0) ? binFall : binFall-1 + (threshold - samplesCorr[binFall-1])/slopeFall;
	trailingEdgeTime = (fittedBinFall + binOffset)*samplingTime;

	// timeOverThreshold
	timeOverThreshold = trailingEdgeTime - leadingEdgeTime;
}

void AhdcExtractor::computeTimeUsingConstantFractionDiscriminator(){
	std::vector<float> signal(binNumber,0.0);
	// signal generation
	for (int bin = 0; bin < binNumber; bin++){
		signal[bin] = (1 - fractionCFD)*samplesCorr[bin]; // we fill it with a fraction of the original signal
		if (bin < binNumber - binDelayCFD)
			signal[bin] += -1*fractionCFD*samplesCorr[bin + binDelayCFD]; // we advance and invert a complementary fraction of the original signal and superimpose it to the previous signal
	}
	// determine the two humps
	int binHumpSup = 0;
	int binHumpInf = 0;
	for (int bin = 0; bin < binNumber; bin++){
		if (signal[bin] > signal[binHumpSup])
			binHumpSup = bin;
	}
	for (int bin = 0; bin < binHumpSup; bin++){ // this loop has been added to be sure : binHumpInf < binHumpSup
		if (signal[bin] < signal[binHumpInf])
			binHumpInf = bin;
	}
	// research for zero
	int binZero = 0;
	for (int bin = binHumpInf; bin <= binHumpSup; bin++){
		if (signal[bin] < 0)
			binZero = bin; // last pass below zero
	} // at this stage : binZero < constantFractionTime/samplingTime <= binZero + 1 // constantFractionTime is determined by assuming a linear fit between binZero and binZero + 1
	float slopeCFD = 0;
	if (binZero + 1 <= binNumber)
		slopeCFD = signal[binZero+1] - signal[binZero];
	float fittedBinZero = (slopeCFD == 0) ? binZero : binZero + (0 - signal[binZero])/slopeCFD;
	constantFractionTime = (fittedBinZero + binOffset)*samplingTime;
	//
	samplesCFD = signal;
}


void AhdcExtractor::Show(const char * filename){
	double tmin = 0;
	double tmax = 6000;
        double threshold = adcOffset + amplitudeFractionCFA*adcMax;

        // Main graph
        double ymax = 0;
        TGraph* gr1 = new TGraph(binNumber);
        for (int i=0;i<binNumber;i++){
                int adc = samplesCorr.at(i) + adcOffset;
                if (ymax < adc) ymax = adc;
                gr1->SetPoint(i,tmin + i*samplingTime,adc);
        }
        // Graph for filling
	int binNumberOVR = 0; // binNumberOverthresholdCFA
	for (int bin=0; bin< binNumber;bin++){
		if ((bin*samplingTime > leadingEdgeTime) && (bin*samplingTime < trailingEdgeTime))
			binNumberOVR++;
	}
        TGraph* gr2 = new TGraph(binNumberOVR+2);
        gr2->SetPoint(0,leadingEdgeTime,threshold);
        gr2->SetPoint(binNumberOVR+1,trailingEdgeTime, threshold);
	int binRiseCFA = (int) leadingEdgeTime/samplingTime;
        for (int i=1;i<=binNumberOVR;i++){
                gr2->SetPoint(i,tmin+samplingTime*(binRiseCFA+i),samplesCorr.at(binRiseCFA+i)+adcOffset);
        }

        // Plot graph
        TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
        gr1->SetTitle("");
        //gr1->SetTitle(TString::Format("%s : p = #Box MeV, #theta = #Box deg, #phi = #Box deg",pid2name[pid].data()));
        gr1->GetXaxis()->SetTitle("Time (ns)");
        gr1->GetXaxis()->SetTitleSize(0.05);
        gr1->GetYaxis()->SetTitle("Charge (adc)");
        gr1->GetYaxis()->SetTitleSize(0.05);
        gr1->GetYaxis()->SetRangeUser(0,ymax+0.05*ymax);
        gr1->SetMarkerColor(kBlack);
        gr1->SetMarkerSize(5);
        gr1->SetLineColor(kBlue);
        gr1->Draw("APL");
        gr2->SetFillColorAlpha(kGreen, 1.0);
        gr2->Draw("F");

        // View decoding
        TLine* line1 = new TLine(leadingEdgeTime,0,leadingEdgeTime,threshold); line1->SetLineWidth(1); line1->SetLineColor(kRed); line1->SetLineStyle(2); line1->Draw(); // leadingEdgeTime
        TLine* line2 = new TLine(0,threshold,leadingEdgeTime,threshold); line2->SetLineWidth(1); line2->SetLineColor(kRed); line2->SetLineStyle(2); line2->Draw(); // leadingEdgeTime
        TLine* line3 = new TLine(trailingEdgeTime,0,trailingEdgeTime,threshold); line3->SetLineWidth(1); line3->SetLineColor(kRed); line3->SetLineStyle(2); line3->Draw(); // trailingEdgeTime
        TArrow* arrow1 = new TArrow(leadingEdgeTime,threshold,trailingEdgeTime,threshold,0.02,"<>"); arrow1->SetLineWidth(1); arrow1->SetLineColor(kRed); arrow1->Draw(); // timeOverThreshold
        TLine* line4 = new TLine(tmin,adcOffset,tmax,adcOffset); line4->SetLineWidth(1); line4->SetLineColor(kRed); line4->SetLineStyle(2); line4->Draw(); // adcOffset
        TLine* line5 = new TLine(0,adcMax+adcOffset,timeMax,adcMax+adcOffset); line5->SetLineWidth(1); line5->SetLineColor(kRed); line5->SetLineStyle(2); line5->Draw(); // adcMax+adcOffset
        TLine* line6 = new TLine(timeMax,0,timeMax,adcMax+adcOffset); line6->SetLineWidth(1); line6->SetLineColor(kRed); line6->SetLineStyle(2); line6->Draw(); // adcMax+adcOffset

        TLatex data;
        data.SetTextSize(0.03);
        data.SetTextAlign(13);
        data.DrawLatexNDC(0.5,0.8,TString::Format("#bf{#bf{leadingEdgeTime} =  %.2lf ns}",leadingEdgeTime).Data());
        data.DrawLatexNDC(0.5,0.8-0.05,TString::Format("#bf{#bf{timeOverThreshold} =  %.2lf ns}",timeOverThreshold).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*2,TString::Format("#bf{#bf{adcMax+adcOffset} =  %.0lf adc }",adcMax+adcOffset).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*3,TString::Format("#bf{#bf{integral} =  %.0lf adc per 44 ns}",integral).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*4,TString::Format("#bf{#bf{adcOffset} =  %d adc}",adcOffset).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*5,"#bf{1 adc =  10^{-5} keV/ns }");
        data.DrawLatexNDC(0.5,0.8-0.05*6,TString::Format("#bf{#bf{constantFractionTime} = %.2lf ns}",constantFractionTime));
        data.SetTextAlign(11);
        data.DrawLatex(leadingEdgeTime,0+(adcMax+adcOffset)*0.02,"leadingEdgeTime");
        data.DrawLatex(timeMax,(adcMax+adcOffset)+(adcMax+adcOffset)*0.02,"adcMax+adcOffset");
        data.DrawLatex(timeMax,threshold,"timeOverThreshold");
        data.DrawLatex(0+tmax*0.02,threshold,"threshold");
        data.DrawLatex(tmax, adcOffset+(adcMax+adcOffset)*0.02,"adcOffset");

        canvas1->Print(filename);
        delete line1; delete line2; delete line3; delete line5;
        delete gr1; delete gr2; delete arrow1;
        delete canvas1;
}

void AhdcExtractor::ShowCFD(const char * filename){
	double tmin = 0;
	double tmax = 6000;
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	TGraph* gr1 = new TGraph(binNumber);
	TGraph* gr2 = new TGraph(binNumber);
	for (int i=0;i<binNumber;i++){
		gr1->SetPoint(i,tmin + i*samplingTime,samplesCorr.at(i)+adcOffset);
		gr2->SetPoint(i,tmin + i*samplingTime,samplesCFD.at(i));
	}
	gr2->SetTitle(TString::Format("#bf{CFD},  fraction = %.1lf, delay = %d index units",fractionCFD,binDelayCFD));
        gr2->GetXaxis()->SetTitle("Time (ns)");
        gr2->GetXaxis()->SetTitleSize(0.05);
        gr2->GetYaxis()->SetTitle("Charge (adc)");
        gr2->GetYaxis()->SetTitleSize(0.05);
        //gr2->SetLineStyle(1);
        gr2->SetLineColor(kRed);
        gr2->SetMarkerColor(kRed);
        gr2->SetMarkerSize(5);
        gr2->GetYaxis()->SetRangeUser(-fractionCFD*adcMax-0.1*adcMax,adcOffset+adcMax+0.1*adcMax);
        gr2->Draw("ALP");

        gr1->SetMarkerColor(kBlue);
        gr1->SetMarkerSize(5);
        gr1->SetLineColor(kBlue);
        //gr1->SetLineStyle(2);
        gr1->Draw("PL");

        TGaxis* axis1 = new TGaxis(tmin,0,tmax,0,tmin,tmax,510,"");
        axis1->SetLineColor(kGreen);
        axis1->SetLabelColor(kGreen);
        //axis1->SetTitle("index units");
        axis1->Draw();

        TLegend* legend = new TLegend(0.7,0.8,0.9,0.9);
        legend->AddEntry(gr1,"Digitized signal","l");
        legend->AddEntry(gr2,"CFD signal","l");
        legend->Draw();

        TLatex data;
        data.SetTextSize(0.04);
        data.SetTextAlign(13);
        data.DrawLatexNDC(0.7,0.6,TString::Format("#bf{#bf{constantFractionTime} =  %.2lf ns}",constantFractionTime).Data());


        //canvas1->Print(TString::Format("./output/SignalCFD_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	canvas1->Print(filename);
        delete gr1; delete gr2; delete canvas1;
        delete axis1; delete legend;
}
