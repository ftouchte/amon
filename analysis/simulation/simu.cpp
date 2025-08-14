/***********************************************
 * Simulation analysis	
 *
 * @author Felix Touchte Codjo
 * @date August 11, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>

#include "reader.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "Math/PdfFuncMathCore.h"
#include "TApplication.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

double breitwigner_pdf(double * x, double * par) {
    return par[0]*ROOT::Math::breitwigner_pdf(x[0], par[1], par[2]);
}

int main(int argc, char const *argv[]) {
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/decoded_22092_20k.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/track/D2/22712/rec_clas_022712.evio.00001.hipo";
    // simu
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/extracted_new_simu.hipo";
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;
    // Histograms
    TH1D* H1_time = new TH1D("leadingEdgeTime", "leadingEdgeTime; time (ns); count", 100, 0, 700); 
    TH1D* H1_timeMax = new TH1D("timeMax", "timeMax; timeMax (ns); count", 100, 200, 900); 
    TH1D* H1_deltaTime = new TH1D("deltaTime", "timeMax - leadingEdgeTime (ns); timeMax - leadingEdgeTime (ns); count", 100, 0, 400); 
    TH1D* H1_sigmaTime = new TH1D("sigmaTime", "(timeMax - leadingEdgeTime)/leadingEdgeTime (ns);(timeMax - leadingEdgeTime)/leadingEdgeTime (ns); count", 100, 0, 1); 
    TH1D* H1_tot = new TH1D("timeOverThreshold", "timeOverThreshol (ns); timeOverThreshol (ns); count", 100, 150, 752); 
    TH1I* H1_wfType = new TH1I("wfType", "wfType; count;", 6, 0, 6); 
    TH1I* H1_adc = new TH1I("amplitude", "amplitude (adc); count;", 100, 0, 3700); 
    TH2D* H2_times = new TH2D("timeMax, leadingEdgeTime", "timeMax vs leadingEdgeTime; timeMax (ns); leadingEdgeTime (ns);", 10, 200, 900, 100, 0, 700); 
    TH2D* H2_amp_tot = new TH2D("amp, tot", "amplitude vs timeOverThreshold;amp (adc); tot (ns)", 100, 0, 3700, 100, 150, 752); 
    TH2D* H2_deltaTime_adc = new TH2D("deltaTime_adc", "deltaTime vs amplitude;deltaTime (ns); amplitude (ns)", 100, 0, 400, 100, 0, 3700); 

    // Loop over events
    while( reader.next()){
        nevents++;
        // Progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(trackBank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);
        
        /*std::vector<int> TrackId;
        for (int i = 0; i < trackBank.getRows(); i++) {
            int trackid = trackBank.getInt("trackid", i);
            if (trackid > 0) {
                TrackId.push_back(trackid);
            }
        }
        if ((int) TrackId.size() < 0) continue;*/
        std::vector<int> HitId;
        for (int i = 0; i < hitBank.getRows(); i++) {
            int trackid = trackBank.getInt("trackid", i);
            int id = hitBank.getInt("id", i);
            if (trackid > 0) {
                HitId.push_back(id);
            }
        }
        if (((int) HitId.size() < 0) && !(mcBank.getRows() > 0)) continue; 
        // Loop over waveforms
        for (int i = 0; i < adcBank.getRows(); i++) {
        //for (int i : HitId) {
            //if (!(adcBank.getInt("wfType", i) <= 1) && !(mcBank.getRows() > 0)) continue;
            double time = adcBank.getFloat("leadingEdgeTime", i);
            double timeMax = adcBank.getFloat("time", i);
            double tot = adcBank.getFloat("timeOverThreshold", i);
            //int wfType = adcBank.getInt("wfType", i);
            int adc = adcBank.getInt("ADC", i);
            //if ((time < 200) || (time > 450)) continue;
            H2_amp_tot->Fill(adc, tot);
            H1_time->Fill(time);
            H1_timeMax->Fill(timeMax);
            H2_times->Fill(timeMax, time);
            H1_deltaTime->Fill(timeMax-time);
            H1_sigmaTime->Fill((timeMax-time)/time);
            H1_tot->Fill(tot);
            //H1_wfType->Fill(wfType);
            H1_adc->Fill(adc);
            H2_deltaTime_adc->Fill(timeMax-time, adc);
        }
        // end Loop over waveforms
    }
    // Process H2_deltaTime_adc
    TCanvas* canvas1 = new TCanvas("c1_deltaTime_amplitude", "deltaTime vs amplitude; deltaTime (ns); amplitude (adc)");
    canvas1->SetLogz();
    H2_deltaTime_adc->Draw("colz");
    int nbins = H2_deltaTime_adc->GetXaxis()->GetNbins();
    printf("> Analyse correlation between deltaTime (ns) and amplitude (adc)\n");
    printf(">   nbins x axis : %d\n", nbins);
    int proj_nbins = 0.05*nbins;
    int Npts = nbins/proj_nbins;
    TGraphErrors* gr = new TGraphErrors(Npts);
    for (int i = 0; i < Npts; i++) {
        TH1D* H1_projTime = H2_deltaTime_adc->ProjectionX(TString::Format("_px_%d", i).Data(), i*proj_nbins, (i+1)*proj_nbins);
        double mean = H1_projTime->GetMean();
        double stdDev = H1_projTime->GetStdDev();
        printf(">   mean : %lf, stdDev : %lf\n", mean, stdDev);
        double y_bin_inf = H2_deltaTime_adc->GetYaxis()->GetBinCenter(i*proj_nbins); // adc_inf for this collection of bins
        double y_bin_sup = H2_deltaTime_adc->GetYaxis()->GetBinCenter((i+1)*proj_nbins); // adc_sup this collection of bins
        gr->SetPoint(i, mean, 0.5*(y_bin_inf+y_bin_sup));
        gr->SetPointError(i, stdDev, 0);
    }
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(21);
    gr->SetLineWidth(2);
    gr->SetLineColor(kRed);
    gr->Draw("same pl");
    TLegend* legend = new TLegend(0.1,0.7,0.25,0.9);
    //legend->SetHeader("deltaTime projection","C");
    legend->AddEntry(gr, "deltaTime projection");
    legend->Draw();

    // end H2_deltaTime_adc analysis 

    // Analysis H1_tot
    /*TCanvas* canvas2 = new TCanvas("c2_tot_fit", "timeOverThreshold (ns); timeOverThreshold (ns); count");
    H1_tot->Draw();
    TF1* func = new TF1("fit_bw_dist", &breitwigner_pdf, 0, 750, 3);
    //TF1* func = new TF1("fit_bw_dist", [&](double*x, double *p){ return p[0]*ROOT::Math::breitwigner_pdf(x[0], p[1], p[2]); }, 0, 750, 3);
    func->SetParameter(0, 18000);
    func->SetParameter(1, 100);
    func->SetParameter(2, 600);
    H1_tot->Fit("fit_bw_dist", "R");*/
    // end H1_tot analysis
    

    printf("nevents    : %ld \n", nevents);
    // output
    //const char * output = "./real_data.root";
    const char * output = "./new_simu_data.root";
    //const char * output = "./real_track_data_after_cuts.root";
    TFile *f = new TFile(output, "RECREATE");
    H1_time->Write("leadingEdgeTime");
    H1_timeMax->Write("timeMax");
    H1_deltaTime->Write("deltaTime");
    H1_sigmaTime->Write("sigmaTime");
    H2_times->Write("timeMax_vs_leadingEdgeTime");
    H1_adc->Write("amplitude");
    H1_tot->Write("timeOverThreshold");
    H2_amp_tot->Write("amplitude_timeOverTheshold");
    H1_wfType->Write("wfType");
    H2_deltaTime_adc->Write("deltaTime_adc");
    canvas1->Write("c1_deltaTime_adc");
    //canvas2->Write("tot_fit");
    f->Close();
    printf("File created : %s\n", output);

    ///////////////////////////////
    // Test
    ///////////////////////////////
    /*TApplication app("ROOT App", nullptr, nullptr);
    TCanvas* canvas3 = new TCanvas("c3_landau", "Landau standard; time; count");
    TF1* f1 = new TF1("f1", [&](double*x, double *p){ return ROOT::Math::landau_pdf(x[0]); }, -5, 15, 0);
    f1->Draw();
    TF1 *f2 = new TF1("f2", "[0]", -5, 15);
    double L_max = f1->GetMaximum(-5, 15);
    f2->SetParameter(0, 0.5*L_max);
    printf("L_max : %lf\n", L_max);
    double t_max = f1->GetMaximumX();
    double t1 = f1->GetX(0.5*L_max,-5, t_max);
    double t2 = f1->GetX(0.5*L_max, t_max, 15);
    printf("t_max : %lf\n", t_max);
    printf("   t1 : %lf\n", t1);
    printf("   t2 : %lf\n", t2);
    printf("t2-t1 : %lf\n", t2-t1);
    printf("t_max-t1 : %lf\n", t_max-t1);
    f2->SetLineColor(kBlue);
    f2->Draw("same l");
    canvas3->Update();
    app.Run();*/
	

    return 0;
}

void progressBar(int state, int bar_length) { // state is a number between 0 and 100
    // for the moment the bar length is not variable
    if (state > bar_length) {return ;}
    printf("\rProgress \033[32m\[");
    for (int i = 0; i <= state; i++) {
        printf("#");
    }
    printf("\033[0m");
    for (int i = state+1; i < bar_length; i++) {
        printf(".");
    }
    if (state == 100) {
        printf("\033[32m] \033[1m %d %%\033[0m\n", state);
    } else {
        printf("] %d %%", state);
    }
    fflush(stdout);
}
