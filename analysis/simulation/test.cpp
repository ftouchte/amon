/***********************************************
 * Test Landau ROOT/CLHEP and CCDB connection	
 *
 * @author Felix Touchte Codjo
 * @date August 16, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
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

#include "CLHEP/GenericFunctions/Landau.hh"

#include <CCDB/Calibration.h>
#include <CCDB/CalibrationGenerator.h>
#include <CCDB/SQLiteCalibration.h>

#include "AhdcCCDB.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

double breitwigner_pdf(double * x, double * par) {
    return par[0]*ROOT::Math::breitwigner_pdf(x[0], par[1], par[2]);
}

int main(int argc, char const *argv[]) {
    TApplication app("ROOT App", nullptr, nullptr);
    /////////////////////////
    /// CCDB
    /// /////////////////////

    TCanvas* canvas2 = new TCanvas("c2_ccdb", "c2_ccdb");
    TH1D* H1_t0 = new TH1D("t0", "t0", 100, 150, 400);
    AhdcCCDB ahdcc;
    for (int wire = 0; wire < 576; wire++) {
        ahdcT0 obj = ahdcc.get_t0(wire);
        H1_t0->Fill(obj.t0);
    }
    H1_t0->Draw();

    ////////////////////////////
    // Landau distribution
    // /////////////////////////
    TCanvas* canvas3 = new TCanvas("c3_landau", "Landau standard");
    TLegend* legend = new TLegend(0.6,0.5,0.95,0.9);
    double sigma = 1;
    std::vector<double> Mu = {0, 1, 2, 4, 8, 11, 15};
    double mu_min = -5*sigma;
    double mu_max = (Mu[Mu.size()-1]+25)*sigma;
    legend->AddEntry("", TString::Format("#sigma = %.2lf", sigma).Data());
    for (int i = 0; i < (int) Mu.size(); i++) {
        //TF1* f1 = new TF1("f1", [&](double*x, double *p){ return ROOT::Math::landau_pdf(x[0], p[0], p[1]); }, mu_min, mu_max, 2);
        TF1* f1 = new TF1("f1", [&](double*x, double *p){ 
            Genfun::Landau L;
            L.peak()  = Genfun::Parameter("Peak", p[1], mu_min, mu_max);
            L.width() = Genfun::Parameter("Width", p[0], 0, 10);
            return L(x[0]);
        }, mu_min, mu_max, 2);
        f1->SetParameter(0, sigma); // stdDev    
        f1->SetParameter(1, Mu[i]);
        f1->SetLineColor(1+i);
        f1->SetTitle("Landau CLHEP");
        //f1->SetTitle("Landau ROOT");
        if (i == 0) { 
            f1->Draw();
        } else {
            f1->Draw("same");
        }
        legend->AddEntry(f1, TString::Format("#mu = %.2lf, x_max = %.3lf", Mu[i], f1->GetMaximumX()).Data());
        printf("%d) Landau x_max : %lf\n", i+1, f1->GetMaximumX());
    }
    legend->Draw();

    printf("///////////////////////////////\n");
    TF1* f1 = new TF1("f1", [&](double*x, double *p){ 
        Genfun::Landau L;
        L.peak()  = Genfun::Parameter("Peak", p[1], -5, 25);
        L.width() = Genfun::Parameter("Width", p[0], 0, 10);
        return L(x[0]);
    }, -5, 25, 2);
    f1->SetParameter(0, 1.0); // stdDev    
    f1->SetParameter(1, 0.0);
    double L_max = f1->GetMaximum(-5, 25);
    double t_max = f1->GetMaximumX();
    double t1 = f1->GetX(0.5*L_max,-5, t_max);
    double t2 = f1->GetX(0.5*L_max, t_max, 15);
    double t10 = f1->GetX(0.1*L_max, -5, t_max);
    double t90 = f1->GetX(0.9*L_max, -5, t_max);
    printf("L_max : %lf\n", L_max);
    printf("t_max : %lf\n", t_max);
    printf("   t1 : %lf\n", t1);
    printf("   t2 : %lf\n", t2);
    printf("t2-t1 : %lf\n", t2-t1);
    printf("tm-t1 : %lf\n", t_max-t1);
    printf("t_10%% : %lf\n", t10);
    printf("t_90%% : %lf\n", t90);
    printf("t_90%%-t_10%% : %lf\n", t90-t10);
    printf("///////////////////////////////\n");

    canvas3->Update();
///////////////////////////////////////
/// Doca 
//////////////////////////////////////
    
    TCanvas* canvas4 = new TCanvas("c4_doca", "Doca from simu");
    TF1* f4 = new TF1("f4", "7*x + 7*x*x + 4*x*x*x", 0, 3);
    f4->Draw();
    canvas4->Update();
    
    TCanvas* canvas5 = new TCanvas("c5_inv_doca", "Doca from simu");
    int Npts4 = 100;
    TGraph* gr4 = new TGraph(Npts4);
    for (int i = 0; i < Npts4; i++) {
        double x = i*3.0/Npts4;
        gr4->SetPoint(i, f4->Eval(x), x);
    }
    gr4->SetLineColor(kBlue);
    gr4->Draw("APL");
    TF1* f5 = new TF1("f5", [&](double*x, double *p){ return p[0] + p[1]*x[0] + p[2]*sqrt(x[0]) + p[3]*pow(x[0], 1.0/3); }, 0, 200, 4);
    gr4->Fit("f5");
    canvas5->Update();


    app.Run();
	

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
