/***********************************************
 * ROOT fitter from data
 *
 * e.g Polynomial fit for time2distance
 *
 * @author Felix Touchte Codjo
 * @date October 30, 2025
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

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

int main(int argc, char const *argv[]) {
    const char * filename = "/home/touchte-codjo/Desktop/amon/analysis/kalman-filter/plot-data-d2t.txt"; 
    std::ifstream flux(filename);
    int Npts = 0;
    if (!flux) {
        std::cout << "Unable to open this file: " << filename << std::endl;
        return 1;
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
    std::cout << "Number of points : " << Npts << std::endl;
    flux.close();
    
    // test
    /*vecx.clear(); vecy.clear();
    for (int i = 0; i < 20; i++) {
        double x = (4.0*i)/20;
        double y = 7*x + 7*x*x + 4*x*x*x;
        vecx.push_back(x);
        vecy.push_back(y);
    }*/



    // xmin, xmax
    double xmin = vecx[0], xmax = vecx[0];
    for (double x : vecx) {
        if (xmin > x) xmin = x;
        if (xmax < x) xmax = x; 
    }
    double xinf = xmin - 0.05*(xmax - xmin);
    double xsup = xmax + 0.05*(xmax - xmin);
    //double xinf = -0.3; 
    //double xsup = 100;
    // Fit data
    TGraph * gr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        gr->SetPoint(i, vecx[i], vecy[i]);
    }
    // template of the fit
    TF1* f1 = new TF1("f1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0) + p[2]*pow(x[0],2.0) + p[3]*pow(x[0],3.0) + p[4]*pow(x[0],4.0) + p[5]*pow(x[0],5.0); }, xinf, xsup, 6); // 6 is the number of parameter
    f1->FixParameter(0,0); // fix the first parameter at 0
    f1->SetParameter(1,1);
    f1->SetParameter(2,1);
    f1->SetParameter(3,1);
    f1->SetParameter(4,1);
    f1->SetParameter(5,1);
    gr->Fit("f1");
    std::vector<double> pars;
    pars.push_back(f1->GetParameter(0));
    pars.push_back(f1->GetParameter(1));
    pars.push_back(f1->GetParameter(2));
    pars.push_back(f1->GetParameter(3));
    pars.push_back(f1->GetParameter(4));
    pars.push_back(f1->GetParameter(5));
    gr->SetTitle(TString::Format("f(x) = %+lf %+lf x %+lf x^{2} %+lf x^{3} %+lf x^{4} %+lf x^{5}; x; y", pars[0],  pars[1], pars[2], pars[3], pars[4], pars[5]).Data());
    
    //////////////////////
    // inverse function
    // ///////////////////
    // ymin, ymax
    double ymin = vecy[0], ymax = vecy[0];
    for (double y : vecy) {
        if (ymin > y) ymin = y;
        if (ymax < y) ymax = y; 
    }
    double yinf = ymin - 0.05*(ymax - ymin);
    double ysup = ymax + 0.05*(ymax - ymin);
    //double yinf = -10; 
    //double ysup = 4;
    TGraph * igr = new TGraph(Npts);
    for (int i = 0; i < Npts; i++) {
        igr->SetPoint(i, vecy[i], vecx[i]);
    }
    // template of the fit
    //TF1* if1 = new TF1("if1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0) + p[2]*pow(x[0],2.0) + p[3]*pow(x[0],3.0) + p[4]*pow(x[0],4.0) + p[5]*pow(x[0],5.0); }, yinf, ysup, 6); // 6 is the number of parameter
    //TF1* if1 = new TF1("if1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0) + p[2]*pow(x[0],1/2.0) + p[3]*pow(x[0],1/3.0) + p[4]*pow(x[0],1/4.0) + p[5]*pow(x[0],1/5.0); }, yinf, ysup, 6); // 6 is the number of parameter
    TF1* if1 = new TF1("if1", [&](double*x, double *p){ return p[0] + p[1]*pow(x[0],1.0) + p[2]*log(1+x[0]) + p[3]*pow(x[0],1/2.0) + p[4]*pow(x[0],1/4.0) + p[5]*pow(x[0],1/5.0); }, yinf, ysup, 6); // 6 is the number of parameter
    if1->FixParameter(0,0); 
    if1->SetParameter(1,0);
    if1->SetParameter(2,0);
    if1->SetParameter(3,0);
    if1->FixParameter(4,0);
    if1->FixParameter(5,0);
    igr->Fit("if1");
    std::vector<double> ipars;
    ipars.push_back(if1->GetParameter(0));
    ipars.push_back(if1->GetParameter(1));
    ipars.push_back(if1->GetParameter(2));
    ipars.push_back(if1->GetParameter(3));
    ipars.push_back(if1->GetParameter(4));
    ipars.push_back(if1->GetParameter(5));
    //igr->SetTitle(TString::Format("f(x) = %+lf %+lf x %+lf x^{2} %+lf x^{3} %+lf x^{4} %+lf x^{5}; x; y", ipars[0],  ipars[1], ipars[2], ipars[3], ipars[4], ipars[5]).Data());
    //igr->SetTitle(TString::Format("f(x) = %+lf %+lf x %+lf x^{1/2} %+lf x^{1/3} %+lf x^{1/4} %+lf x^{1/5}; x; y", ipars[0],  ipars[1], ipars[2], ipars[3], ipars[4], ipars[5]).Data());
    igr->SetTitle(TString::Format("f(x) = %+lf %+lf x %+lf ln(1+x) %+lf x^{1/2} %+lf x^{1/4} %+lf x^{1/5}; x; y", ipars[0],  ipars[1], ipars[2], ipars[3], ipars[4], ipars[5]).Data());

    // test
    TF1* f2 = new TF1("f2", "7*pow(x,1.0) + 7*pow(x,2.0) + 4*pow(x,3.0)", 0, 4);

    const char * output = "./output/fitter.root";
    TFile *f = new TFile(output, "RECREATE");
    gr->Write("graph");
    igr->Write("inv_graph");
    f2->Write("test_d2t_old");
    //canvas->Write("canvas");
    f->Close();
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

