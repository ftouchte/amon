/***********************************************
 * Study the Kalman Filter for the AHDC
 *
 * Distance in simulation
 *
 * @author Felix Touchte Codjo
 * @date October 24, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

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

int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
    //const char * output = "./output/distance.root";
    const char * output = "./output/distance-v15.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v15.hipo";
    //const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v7.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    //hipo::bank  trackBank(factory.getSchema("AHDC::track"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nMCtracks =0;
    long unsigned int nKFtracks =0;
    // Histograms
    // hit 
    TH1D* H1_mcDocaTime  = new TH1D("mcDocaTime", "mcDocaTime; mcDocaTime (ns); #count", 100, 0, 250); 
    TH1D* H1_decodedTime = new TH1D("decodedTime", "decodedTime; decodedTime (ns); count", 100, 0, 250);
    TH1D* H1_deltaTime = new TH1D("deltaTime", "deltaTime = mcDocaTime - decodedTime; deltaTime (ns); count", 100, -50, 50);
    TH1D* H1_mcDoca    = new TH1D("mcDoca", "mcDoca; mcDoca (mm); #count", 100, 0, 4); 
    TH1D* H1_decodedDistance = new TH1D("decodedDistance", "decodedDistance; decodedDistance (mm); count", 100, 0, 4);
    TH1D* H1_deltaDistance = new TH1D("deltaDistance", "deltaDistance = mcDoca - decodedDistance; deltaDistance (mm); count", 100, -1, 1);
    TH1D* H1_mcMeanTime  = new TH1D("mcMeanTime", "mcMeanTime; mcMeanTime (ns); count", 100, 0, 250);
    TH1D* H1_mcMeanDistance  = new TH1D("mcMeanDistance", "mcMeanDistance; mcMeanDistance (mm); count", 100, 0, 4);
    std::vector<TH2D*> H2_corrDistances;
    H2_corrDistances.push_back(new TH2D("corr_1", "mcDoca vs decodedDistance;mcDoca (mm); decodedDistance (mm)", 100, 0, 4, 100, 0, 4));
    H2_corrDistances.push_back(new TH2D("corr_2", "mcDoca vs mcMeanDistance;mcDoca (mm); mcMeanDistance (mm)", 100, 0, 4, 100, 0, 4));
    H2_corrDistances.push_back(new TH2D("corr_3", "mcMeanDistance vs decodedDistance;mcMeanDistance (mm); decodedDistance (mm)", 100, 0, 4, 100, 0, 4));
    std::vector<TH2D*> H2_corrTimes;
    H2_corrTimes.push_back(new TH2D("corr_4", "mcDocaTime vs decodedTime;mcDocaTime (ns); decodedTime (ns)", 100, 0, 250, 100, 0, 250));
    H2_corrTimes.push_back(new TH2D("corr_5", "mcDocaTime vs mcMeanTime;mcDocaTime (ns); mcMeanTime (ns)", 100, 0, 250, 100, 0, 250));
    H2_corrTimes.push_back(new TH2D("corr_6", "mcMeanTime vs decodedTime;mcMeanTime (ns); decodedTime (ns)", 100, 0, 250, 100, 0, 250));
    //TH2D* H2_time2distance = new TH2D("corr_t2d", "mcMeanTime vs mcDoca;mcMeanTime (ns); mcDoca (mm)", 100, 0, 250, 100, 0, 4);
    //TH1D* H1_t0  = new TH1D("t0", "t0; t0 (ns); count", 100, 0, 400);
    //TH1D* H1_remcDoca    = new TH1D("remcDoca", "remcDoca from t2d; remcDoca (mm); #count", 100, 0, 4); 
    TH1D* H1_adc    = new TH1D("adc", "ADC; ADC (adc); #count", 100, 0, 4000); 
    TH1D* H1_tot    = new TH1D("tot", "timeOverThreshold; timeOverThreshold (ns); #count", 100, 150, 750); 
    TH2D* H2_deltaTime_adc = new TH2D("corr_deltaTime_adc", "deltaTime vs ADC;ADC (adc); deltaTime (ns)", 100, 0, 4000, 100, -50, 50);
    TH2D* H2_deltaDistance_adc = new TH2D("corr_deltaDistance_adc", "deltaDistance vs ADC;ADC (adc); deltaDistance (mm)", 100, 0, 4000, 100, -2.6, 2.6);
    AhdcCCDB ahdcConstants;
    //****** T2D inversion************************
    ahdcT2d T2D = ahdcConstants.get_t2d();
    auto eval_t2d = [T2D] (double x) -> double {
        return T2D.p0 + T2D.p1*pow(x, 1.0) + T2D.p2*pow(x, 2.0) + T2D.p3*pow(x, 3.0) + T2D.p4*pow(x, 4.0) + T2D.p5*pow(x, 5.0);
    };
    double xi[50];
    double yi[50];
    for (int i = 0; i < 50; i++) {
        xi[i] = i*(350.0/50);
        yi[i] = eval_t2d(xi[i]);
    }
    // inverse of the xime2yistance
    auto eval_inv_t2d = [xi, yi] (double y) -> double {
        if (y < 0) {
            return ((xi[1]-xi[0])/(yi[1]-yi[0]))*(y - yi[0]) + xi[0];
        }
        else if (y >= yi[49]) {
            return ((xi[49]-xi[48])/(yi[49]-yi[48]))*(y - yi[48]) + xi[48];
        } else {
            int i = 0;
            while (i < 48) {
                if ((y >= yi[i]) && (y < yi[i+1])) {
                        break;
                }
                i++;
            }
            return ((xi[i+1]-xi[i])/(yi[i+1]-yi[i]))*(y - yi[i]) + xi[i];
        }
    };
    //*********** END T2D inversion ****************
    int Npts = 100;
    TGraph * gr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        double x = i*(350.0/Npts);
        double y = eval_t2d(x);
        gr->SetPoint(i, x, y);
    }
    TGraph * igr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        double x = i*(4.0/Npts);
        double y = eval_inv_t2d(x);
        igr->SetPoint(i, x, y);
    }
    gr->SetTitle("time2distance; time (ns); distance (mm)");
    igr->SetTitle("distance2time; distance (mm); time (ns)");
    TGraph * igr_gr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        double x = i*(350.0/Npts);
        double y = eval_t2d(x);
        double xx = eval_inv_t2d(y);
        igr_gr->SetPoint(i, x, xx);
    }
    TGraph * gr_igr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        double x = i*(4.0/Npts);
        double y = eval_inv_t2d(x);
        double xx = eval_t2d(y);
        gr_igr->SetPoint(i, x, xx);
    }
    igr_gr->SetTitle("igr o gr; time (ns); time (ns)");
    gr_igr->SetTitle("gr o igr; distance (mm); distance (mm)");


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

        nMCtracks += mcBank.getRows();
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
                double mcMeanTime = 0.01*wfBank.getInt(std::string("s30").c_str(), i);
                double mcDoca = 0.001*wfBank.getInt(std::string("s28").c_str(), i);
                //double p0 = 0;
                //double p1 = 2.55132;
                //double p2 = 10.7884;
                //double p3 = 12.8042;
                //double p4 = -9.91149;
                //double p5 = 2.38082;
                //double mcDocaTime = p0 + p1*mcDoca + p2*pow(mcDoca,2) + p3*pow(mcDoca,3) + p4*pow(mcDoca,4) + p5*pow(mcDoca,5); 
                double mcDocaTime = eval_inv_t2d(mcDoca);
                //double t0 = ahdcConstants.get_t0(sector, layer, component).t0;
                //double decodedTime = adcBank.getDouble("leadingEdgeTime", i) - t0; // for the simu, startTime = 0
                double decodedTime = hitBank.getDouble("time", h); 
                double decodedDistance = hitBank.getDouble("doca", h);
                //double mcMeanDistance = (mcMeanTime < 0) ? 0 : -0.0497 -0.00667*mcMeanTime + 0.389*sqrt(mcMeanTime) - 0.189*pow(mcMeanTime, 1.0/3);
                double mcMeanDistance = -99;
                H1_mcMeanTime->Fill(mcMeanTime);
                H1_decodedTime->Fill(decodedTime);
                H1_mcDocaTime->Fill(mcDocaTime);
                H1_deltaTime->Fill(mcDocaTime - decodedTime);
                H1_deltaDistance->Fill(mcDoca - decodedDistance);
                H1_mcMeanDistance->Fill(mcMeanDistance);
                H1_decodedDistance->Fill(decodedDistance);
                H1_mcDoca->Fill(mcDoca);
                //H1_t0->Fill(t0);
                H2_corrDistances[0]->Fill(mcDoca, decodedDistance);
                H2_corrDistances[1]->Fill(mcDoca, mcMeanDistance);
                H2_corrDistances[2]->Fill(mcMeanDistance, decodedDistance);
                H2_corrTimes[0]->Fill(mcDocaTime, decodedTime);
                H2_corrTimes[1]->Fill(mcDocaTime, mcMeanTime);
                H2_corrTimes[2]->Fill(mcMeanTime, decodedTime);
                //H2_time2distance->Fill(mcMeanTime, mcDoca);
                //H1_remcDoca->Fill(-0.0129887*mcDocaTime -0.331228*log(1+mcDocaTime) +0.516754*sqrt(mcDocaTime));
                int adc = adcBank.get("ADC",i);
                H1_adc->Fill(adc);
                H2_deltaTime_adc->Fill(adc, mcDocaTime - decodedTime);
                H2_deltaDistance_adc->Fill(adc, mcDoca - decodedDistance);
                H1_tot->Fill(adcBank.get("timeOverThreshold",i));
            }
        }
    }
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    TFile *f = new TFile(output, "RECREATE");
    
    //H1_mcMeanTime->Write("mcMeanTime");
    H1_decodedTime->Write("decodedTime");
    H1_mcDocaTime->Write("mcDocaTime");
    //H1_mcMeanDistance->Write("mcMeanDistance");
    H1_decodedDistance->Write("decodedDistance");
    H1_mcDoca->Write("mcDoca");
    H2_corrDistances[0]->Write("corr_mcDoca_decodedDistance");
    //H2_corrDistances[1]->Write("corr_mcDoca_mcMeanDistance");
    //H2_corrDistances[2]->Write("corr_mcMeanDistance_decodedDistance");
    H2_corrTimes[0]->Write("corr_mcDocaTime_decodedTime");
    //H2_corrTimes[1]->Write("corr_mcDocaTime_mcMeanTime");
    //H2_corrTimes[2]->Write("corr_mcMeanTime_decodedTime");
    //H2_time2distance->Write("time2distance");
    //H1_t0->Write("t0");
    H1_deltaTime->Write("deltaTime");
    H1_deltaDistance->Write("deltaDistance");
    //H1_remcDoca->Write("remcDoca");
    H1_adc->Write("ADC");
    H2_deltaTime_adc->Write("corr_deltaTime_adc");
    H2_deltaDistance_adc->Write("corr_deltaDistance_adc");
    H1_tot->Write("ToT");
    gr->Write("gr");
    igr->Write("igr");
    igr_gr->Write("igr_gr");
    gr_igr->Write("gr_igr");
    TCanvas *c11 = new TCanvas("c11", "time2distance funtions");
    c11->Divide(2,2);
    c11->cd(1);
    gr->Draw("APL");
    c11->cd(2);
    igr->Draw("APL");
    c11->cd(3);
    igr_gr->Draw("APL");
    TF1 * f1 = new TF1("f1", "x", 0, 350);
    f1->Draw("same");
    c11->cd(4);
    gr_igr->Draw("APL");
    TF1 * f2 = new TF1("f2", "x", 0, 4);
    f2->Draw("same");
    c11->Write("t2d_functions");

    f->Close();
    printf("File created : %s\n", output);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
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
