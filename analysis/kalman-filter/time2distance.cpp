/***********************************************
 * Study the Kalman Filter for the AHDC
 *
 * time2distance from AHDC::hits
 *
 * @author Felix Touchte Codjo
 * @date October 24, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
//#include <algorithm>

#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

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
#include "TText.h"
#include "THStack.h"

#include "AhdcCCDB.h"
#include "futils.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities
TCanvas* time2distance(std::vector<double> & pars);

int main(int argc, char const *argv[]) {

    // output
    const char * output = "./output/time2distance-v10-bis.root";
    TFile *f = new TFile(output, "RECREATE");
    TDirectory *projection_dir = f->mkdir("t2d_fit_projection1D");
    //std::vector<double> pars;
    //TCanvas * canvas1 = time2distance(pars);
    //if (canvas1) canvas1->Write();
    //f->Close();
    //return 0;
    auto start = std::chrono::high_resolution_clock::now();
    
    
    /////////////////////////
    /// simu
    /// /////////////////////
    //const char* filename = "/home/touchte-codjo/Desktop/amon/analysis/simulation/elastics_events.hipo";
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/elastics_events-v10.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nKFtracks =0;
    // Histograms
    // hit 
    TH1D* H1_time = new TH1D("time", "time; time (ns); count", 100, 0, 250);
    TH1D* H1_distance = new TH1D("distance", "distance; distance (mm); count", 100, 0, 4);
    TH2D* H2_time2distance = new TH2D("corr_t2d", "time2distance;time (ns); distance (mm)", 100, 0, 320, 100, 0, 4);
    TH1D* H1_residual = new TH1D("residuals", "residuals; residuals (mm); #count", 100, -2, 2);
    TH2D* H2_time2distance_filtered = new TH2D("corr_t2d_filtered", "time2distance;time (ns); distance (mm)", 100, 0, 320, 100, 0, 4);
    TH1D* H1_residual_filtered = new TH1D("residuals_filtered", "residuals; residuals (mm); #count", 100, -2, 2);
    //TH1D* H1_t0  = new TH1D("t0", "t0; t0 (ns); count", 100, 0, 400);
    //TH1D* H1_testDistance = new TH1D("testDistance", "testDistance; testDistance (mm); count", 100, 0, 4);
    //TH1D* H1_leadingEdgeTime = new TH1D("leadingEdgeTime", "leadingEdgeTime; leadingEdgeTime (ns); count", 100, 250, 700);

    AhdcCCDB ahdcConstants;
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

        nKFtracks += trackBank.getRows();
        
        for (int t = 0; t < trackBank.getRows(); t++) {
            for (int h = 0; h < hitBank.getRows(); h++) {
                if (hitBank.getInt("trackid", h) != trackBank.getInt("trackid", t)) { 
                    continue;
                }
                int i = hitBank.getInt("id", h) - 1;
                //int sector = wfBank.getShort("sector", i);
                //int layer = wfBank.getShort("layer", i);
                //int component = wfBank.getShort("component", i);
                //double t0 = ahdcConstants.get_t0(sector, layer, component).t0;
                double time = hitBank.getDouble("time", h); 
                double residual = hitBank.getDouble("residual", h);
                double doca = hitBank.getDouble("doca", h);
                double distance = doca - residual; 
                //double distance = (time < 0) ? 0 : -0.0497 -0.00667*time + 0.389*sqrt(time) - 0.189*pow(time, 1.0/3);
                H1_time->Fill(time);
                H1_distance->Fill(distance);
                //H1_t0->Fill(t0);
                H2_time2distance->Fill(time, distance);
                H1_residual->Fill(residual);
                if (std::abs(residual) > 1e-10) { // remove the peak at 0
                    H2_time2distance_filtered->Fill(time, distance);
                    H1_residual_filtered->Fill(residual);
                }
                //double value = pars[0] + pars[1]*time + pars[2]*pow(time, 1.0/2) + pars[3]*pow(time, 1.0/3);
                //H1_testDistance->Fill(value);
                //H1_leadingEdgeTime->Fill(adcBank.getFloat("leadingEdgeTime", i));
            }
        }
    }
    ///////////////////////////////////////////////////
    // select maximum distance for the time interval
    // ////////////////////////////////////////////////
    TCanvas* canvas2 = new TCanvas("c2_time2distance_fitted", "c2_time2distance_fitted");
    //canvas2->SetLogz();
    H2_time2distance_filtered->Draw("colz");
    TText* text1 = new TText();
    text1->SetTextSize(0.02);
    text1->SetTextColor(kGreen);
    int nbins = H2_time2distance_filtered->GetXaxis()->GetNbins();
    printf("> Fitting of the time2distance\n");
    printf(">   nbins x axis : %d\n", nbins);
    //int proj_nbins = 0.05*nbins;
    int proj_nbins = 2;
    int Npts = nbins/proj_nbins;
    TGraphErrors* gr = new TGraphErrors(Npts);
    for (int i = 0; i < Npts; i++) {
        TCanvas* canvas3 = new TCanvas("c3", "c3_proj");
        canvas3->cd();
        TH1D* H1_projDistance = H2_time2distance_filtered->ProjectionY(TString::Format("_py_%d", i).Data(), i*proj_nbins, (i+1)*proj_nbins);
        H1_projDistance->Fit("gaus");
        projection_dir->cd();
        H1_projDistance->Write(TString::Format("_py_%d", i).Data());
        TF1 *fgaus = H1_projDistance->GetFunction("gaus");
        double mean = fgaus->GetParameter("Mean");
        double stdDev = fgaus->GetParameter("Sigma");
        //double mean = H1_projDistance->GetMean();
        //double stdDev = H1_projDistance->GetStdDev();
        printf(">   mean : %lf, stdDev : %lf\n", mean, stdDev);
        double x_bin_inf = H2_time2distance_filtered->GetXaxis()->GetBinCenter(i*proj_nbins);
        double x_bin_sup = H2_time2distance_filtered->GetXaxis()->GetBinCenter((i+1)*proj_nbins);
        gr->SetPoint(i, 0.5*(x_bin_sup+x_bin_inf), mean);
        gr->SetPointError(i, 0, stdDev);
        canvas2->cd();
        //text1->DrawText(0.5*(x_bin_sup+x_bin_inf), mean, TString::Format("%.2lf +/- %.2lf", mean, stdDev).Data());
    }
    canvas2->cd();
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(21);
    gr->SetLineWidth(1);
    gr->SetLineColor(kGreen);
    gr->Draw("same pl");
    TLegend* legend = new TLegend(0.1,0.7,0.6,0.9);
    //legend->SetHeader("deltaTime projection","C");
    legend->AddEntry(gr, "most probable distance for a given time slot");
    ///////////////////////
    // Fit the points
    // ////////////////////
    TF1* f1 = new TF1("f1", [&](double*x, double *p){ return p[0]*pow(x[0],1.0) + p[1]*pow(x[0],2.0) + p[2]*pow(x[0],3.0); }, 0, 100, 3); // 3 is the number of parameter
    //f1->SetLineColor(kGreen);
    f1->SetLineWidth(2);
    f1->SetParameter(0,1);
    f1->SetParameter(1,1);
    f1->SetParameter(2,1);
    //gr->SetPoint(0, 0, 0); // force the first point to be at zero
    gr->Fit("f1");
    std::vector<double> pars;
    pars.push_back(f1->GetParameter(0));
    pars.push_back(f1->GetParameter(1));
    pars.push_back(f1->GetParameter(2));
    legend->AddEntry(f1, TString::Format("polynomial fit: %+lf#bf{t}%+lf#bf{t}^{2}%+lf#bf{t}^{3}", pars[0], pars[1], pars[2]).Data());
    canvas2->Update();
    legend->Draw();

    f->cd();
    H1_time->Write("time");
    H1_distance->Write("distance");
    H2_time2distance->Write("time2distance");
    H1_residual->Write("residual");
    H2_time2distance_filtered->Write("time2distance_filtered");
    H1_residual_filtered->Write("residual_filtered");
    canvas2->Write("time2distance_fitted");
    //H1_testDistance->Write("testDistance");
    //H1_leadingEdgeTime->Write("leadingEdgeTime");
    //H1_t0->Write("t0");

    f->Close();
    printf("File created : %s\n", output);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("nevents    : %ld \n", nevents);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("* time elapsed : %lf seconds\n", elapsed.count());
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

TCanvas* time2distance(std::vector<double> & pars) {
    // Rea data
    //const char * filename = "/home/touchte-codjo/Desktop/amon/analysis/kalman-filter/plot-data.txt"; 
    const char * filename = "/home/touchte-codjo/Desktop/amon/analysis/kalman-filter/plot-data-2.txt"; 
    std::ifstream flux(filename);
    int Npts = 0;
    if (!flux) {
        std::cout << "Unable to open this file: " << filename << std::endl;
        return nullptr;
    }
    std::string title1, title2;
    flux >> title1;
    flux >> title2;
    std::cout << title1 << ", " << title2 << std::endl;
    double x, y;
    std::vector<double> vecx, vecy;
    while (!flux.eof()) {
        flux >> x;
        flux >> y;
        std::cout << x << ", " << y << std::endl;
        Npts++;
        vecx.push_back(x);
        vecy.push_back(y);
    }
    std::cout << "Number of poiints : " << Npts << std::endl;
    flux.close();
    // Fit data
    TCanvas* canvas = new TCanvas("c1", "c1");
    TGraph * gr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        gr->SetPoint(i, vecx[i], vecy[i]);
    }
    //gr->SetTitle("time2distance (Page 67, Lucien's thesis); time (ns); distance (mm)");
    gr->SetTitle("fit from data; time (ns); distance (mm)");
    gr->Draw("APL");
    //TF1* f1 = new TF1("f1", [&](double*x, double *p){ return p[0]*(p[1] - 1/(x[0] + p[2]*pow(x[0],3) + p[3])); }, 0, 100, 4); // 4 is the number of parameter
    //TF1* f1 = new TF1("f1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0/1) + p[2]*pow(x[0],1.0/2) + p[3]*pow(x[0],1.0/3); }, 0, 100, 4); // 4 is the number of parameter
    TF1* f1 = new TF1("f1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0) + p[2]*pow(x[0],2.0) + p[3]*pow(x[0],3.0); }, 0, 100, 4); // 4 is the number of parameter
    f1->SetParameter(0,1);
    f1->SetParameter(1,1);
    f1->SetParameter(2,1);
    f1->SetParameter(3,1);
    gr->Fit("f1");
    canvas->Update();
    pars.clear();
    pars.push_back(f1->GetParameter(0));
    pars.push_back(f1->GetParameter(1));
    pars.push_back(f1->GetParameter(2));
    pars.push_back(f1->GetParameter(3));
    return canvas;
}







