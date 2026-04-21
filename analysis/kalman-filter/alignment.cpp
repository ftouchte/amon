/**
 * @brief Make summary plots related to the AHDC alignment
 * 
 * Code resuming the AHDC position scan
 * 
 * @date March 31, 2026
 */

#include "alignment.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "Math/PdfFuncMathCore.h"
#include "TText.h"
#include "THStack.h"
#include "TArrow.h"
#include "TLine.h"


int main(int argc, char const *argv[]) {
    TGraphErrors* gre = new TGraphErrors();
    TGraph* g_mean = new TGraph();
    TGraph* g_dev = new TGraph();

    PositionScan scans;
    alignmentEntry optimum({1e10, 1e10, 1e10});
    double min_dev = 1e10;
    for (alignmentEntry entry : scans.values) {
        double clas_alignment = entry.clas_alignment;
        double mean_angle = entry.mean_angle;
        double deviation_angle = entry.deviation_angle;
        gre->AddPointError(-clas_alignment, mean_angle, 0, deviation_angle);
        g_mean->AddPoint(-clas_alignment, mean_angle);
        g_dev->AddPoint(-clas_alignment, deviation_angle);

        if (deviation_angle < min_dev) {
            optimum = entry;
            min_dev = deviation_angle;
        }
    }

    // Title
    gre->SetTitle("Angle dispersion versus AHDC position");
    gre->GetXaxis()->SetTitle("shift (mm)");
    gre->GetYaxis()->SetTitle("angle dispersion (deg)");
    gre->SetMarkerStyle(20);

    g_mean->SetTitle("Mean angle versus AHDC position");
    g_mean->GetXaxis()->SetTitle("shift (mm)");
    g_mean->GetYaxis()->SetTitle("mean angle (deg)");

    g_dev->SetTitle("Angle dispersion versus AHDC position");
    g_dev->GetXaxis()->SetTitle("shift (mm)");
    g_dev->GetYaxis()->SetTitle("angle dispersion (deg)");

    // Main plot

    TGraphErrors* gre_optimum = new TGraphErrors();
    gre_optimum->AddPointError(-optimum.clas_alignment, optimum.mean_angle, 0, optimum.deviation_angle);
    gre_optimum->SetMarkerColor(kRed);
    gre_optimum->SetLineColor(kRed);
    gre_optimum->SetMarkerStyle(20);

    // Legend
    TLegend * legend = new TLegend(0.1,0.7, 0.65, 0.9);
    //legend->AddEntry(gre, "mean angle and dispersion", "lfp");
    legend->AddEntry(gre, "mean angle", "p");
    legend->AddEntry(gre, "angle dispersion", "e");
    legend->AddEntry(gre_optimum, TString::Format("optimal dispersion : (%.2lf mm, %.2lf #pm %.2lf deg)", -optimum.clas_alignment, optimum.mean_angle, optimum.deviation_angle).Data(), "lfp");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    printf("optimum ===> clas alignment : %lf , angle deviation : %lf, mean angle : %lf\n", -optimum.clas_alignment, optimum.mean_angle, optimum.deviation_angle);

    TCanvas* c = new TCanvas("angle_dispersion_vs-position");
    gre->Draw("ap");
    gre_optimum->Draw("same p");
    legend->Draw();

    // Output
    TFile *f = new TFile("./output/alignment.root", "RECREATE");
    //gre->Write("combined_angle_dipersion");
    c->Write("combined_angle_dipersion");
    g_mean->Write("mean_angle");
    g_dev->Write("angle_dispersion");
    f->Close();

    printf("File created : %s\n", "./output/alignment.root");
}

