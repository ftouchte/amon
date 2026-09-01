#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"


void wfType_occupancy_study() {

    TFile *f1 = new TFile("./output/wfType-study-v1.root", "READ");
    TH1F *h1 = (TH1F*)f1->Get("occupancy_from_selected_hits");
    
    TH1F *h0 = (TH1F*)f1->Get("occupancy_from_all_hits");

    TFile *f2 = new TFile("./output/wfType-study-v2.root", "READ");
    TH1F *h2 = (TH1F*)f2->Get("occupancy_from_selected_hits");

    TFile *f3 = new TFile("./output/wfType-study-v3.root", "READ");
    TH1F *h3 = (TH1F*)f3->Get("occupancy_from_selected_hits");

    TFile *f4 = new TFile("./output/wfType-study-v4.root", "READ");
    TH1F *h4 = (TH1F*)f4->Get("occupancy_from_selected_hits");

    TFile *f5 = new TFile("./output/wfType-study-v5.root", "READ");
    TH1F *h5 = (TH1F*)f5->Get("occupancy_from_selected_hits");

    TCanvas* c1 = new TCanvas("c1","c1", 1200, 800);

    h0->SetLineColor(kBlack);
    h1->SetLineColor(kBlack);
    h1->SetLineStyle(2);
    h2->SetLineColor(kGreen);
    h3->SetLineColor(kBlue);
    h4->SetLineColor(kBlack);
    h4->SetLineStyle(7);
    h5->SetLineColor(kRed);
    

    h0->GetXaxis()->SetRange(0,47);
    h0->SetStats(0);
    h0->SetXTitle("wire number");
    h0->SetYTitle("occupancy [%]");
    h0->SetTitle("Occupancy on the first layer");

    h0->Draw("hist");
    h1->Draw("hist same");
    h2->Draw("hist same");
    h3->Draw("hist same");
    h4->Draw("hist same");
    h5->Draw("hist same");

    h0->SetLineWidth(2);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);
    h4->SetLineWidth(2);
    h5->SetLineWidth(2);

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(h0, "All hits",  "lfp");
    legend->AddEntry(h1, "Raw cuts & ADC #geq 200", "lfp");
    legend->AddEntry(h2, "Raw cuts & ADC #geq 50", "lfp");
    legend->AddEntry(h3, "Raw cuts & ADC #geq 0", "lfp");
    legend->AddEntry(h4, "wfType #leq 2 only", "lfp");
    legend->AddEntry(h5, "wfType #leq 2 & Raw cuts & ADC #geq 0", "lfp");
    legend->Draw("same");


    TFile *f = new TFile("./output/wfType_occupancy_study.root", "RECREATE");
    c1->Write("combined_occupancy");
    f->Close();
    c1->Show();

}