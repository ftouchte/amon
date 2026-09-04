#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"


void wfType_occupancy_study() {

    TFile *f1 = new TFile("./output/occupancy-v1.root", "READ");
    TH1F *h1 = (TH1F*)f1->Get("occupancy_from_selected_hits");
    
    TH1F *h0 = (TH1F*)f1->Get("occupancy_from_all_hits");

    TFile *f2 = new TFile("./output/occupancy-v2.root", "READ");
    TH1F *h2 = (TH1F*)f2->Get("occupancy_from_selected_hits");

    TFile *f3 = new TFile("./output/occupancy-v3.root", "READ");
    TH1F *h3 = (TH1F*)f3->Get("occupancy_from_selected_hits");

    TFile *f4 = new TFile("./output/occupancy-v4.root", "READ");
    TH1F *h4 = (TH1F*)f4->Get("occupancy_from_selected_hits");

    TFile *f5 = new TFile("./output/occupancy-v5.root", "READ");
    TH1F *h5 = (TH1F*)f5->Get("occupancy_from_selected_hits");

    TCanvas* c1 = new TCanvas("c1","c1", 1200, 800);

    h0->SetLineColor(kBlack);
    h1->SetLineColor(kBlack);
    h1->SetLineStyle(2);
    h2->SetLineColor(kGreen);
    h3->SetLineColor(kBlue);
    h4->SetLineColor(kBlack);
    h4->SetLineStyle(9);
    h5->SetLineColor(kRed);
    

    h0->GetXaxis()->SetRange(0,47);
    h0->SetStats(0);
    h0->SetXTitle("wire number");
    h0->SetYTitle("occupancy [%]");
    //h0->SetTitle("Occupancy on the first layer");
    h0->SetTitle("");

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

    TLegend* legend = new TLegend(0.55,0.65,0.89,0.89);
    legend->AddEntry(h0, "All hits",  "l");
    legend->AddEntry(h1, "Raw cuts & ADC #geq 200", "l");
    legend->AddEntry(h2, "Raw cuts & ADC #geq 50", "l");
    legend->AddEntry(h3, "Raw cuts & ADC #geq 0", "l");
    legend->AddEntry(h4, "wfType #leq 2 only", "l");
    legend->AddEntry(h5, "wfType #leq 2 & Raw cuts & ADC #geq 0", "l");
    legend->Draw("same");

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h0->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h0 : %lf %%\n", sum);
    }

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h1->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h1 : %lf %%\n", sum);
    }

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h2->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h2 : %lf %%\n", sum);
    }

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h3->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h3 : %lf %%\n", sum);
    }

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h4->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h4 : %lf %%\n", sum);
    }

    {
        double sum = 0;
        for (int i = 1; i <= 47; i++) {
            sum += h5->GetBinContent(i);
        }
        sum /= 47;
        printf("occ h5 : %lf %%\n", sum);
    }


    TFile *f = new TFile("./output/wfType_occupancy_study.root", "RECREATE");
    c1->Write("combined_occupancy");
    f->Close();
    c1->Show();

}