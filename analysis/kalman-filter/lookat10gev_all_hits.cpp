/**********************************************
 * Copy of lookat10gev.cpp
 * 
 * Study time reconstruction of saturated ADC
 * 
 * Here we look at all hits
 * 
 * @author Felix Touchte Codjo
 * @date August 19, 2026
 **********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

#include "lookat10gev_all_hits.h"
#include "fOptions.h"
#include "Units.h"

#include "reader.h"

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
#include "TTree.h"
#include "TMath.h"
#include "TMultiGraph.h"


double ahdcPulse(double x, double mpv, double sigma, double norm, double constant) {
    return constant + norm*TMath::Landau(x, mpv, sigma, true);
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return ahdcPulse(x[0], par[0], par[1], par[2], par[3]);
}

TGraph* getSlopes(TGraph * gr) {
    int N = gr->GetN();
    TGraph * dgr = new TGraph(N-1);
    for (int i = 0; i < N-1; i++) {
        double x1 = gr->GetPointX(i);
        double y1 = gr->GetPointY(i);
        double x2 = gr->GetPointX(i+1);
        double y2 = gr->GetPointY(i+1);
        double slope = (y2-y1)/(x2-x1);
        dgr->SetPoint(i,x1, slope);
    }
    dgr->SetMarkerColor(kRed);
    dgr->SetMarkerStyle(20);
    dgr->SetTitle("Derivate of the pulse; time (ns); slope (ADC/ns)");
    
    return dgr;
}

void fitGaus(TH1D* h) {
    double mean = h->GetMean();
    double sigma = h->GetStdDev();
    h->Fit("gaus", "", "R", mean-sigma, mean+sigma);
}


int main(int argc, char const *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();
    
    /// --- Load files from options
    fOptions OPT({"-i", "-o"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");


    /// --- Histograms
    Histograms* Histos = new Histograms();

    TTree *tree = new TTree("T", "ADC versus Time");

    double residual;
    double time;
    int adc;
    int raw_adc;
    double tot;
    double leadingEdgeTime;
    int wfType;
    double slope_max;
    double time_slope_max;

    tree->Branch("residual", &residual, "residual/D");
    tree->Branch("time", &time, "time/D");
    tree->Branch("leadingEdgeTime", &leadingEdgeTime, "leadingEdgeTime/D");
    tree->Branch("tot", &tot, "tot/D");
    tree->Branch("adc", &adc, "adc/I");
    tree->Branch("raw_adc", &raw_adc, "raw_adc/I");
    tree->Branch("wfType", &wfType, "wfType/I");
    tree->Branch("slope_max", &slope_max, "slope_max/D");
    tree->Branch("time_slope_max", &time_slope_max, "time_slope_max/D");


    int nb_wf0 = 0;
    int nb_wf0_outlier = 0;
    int nb_wf0_truncated = 0;
    int nb_wf0_functional = 0;
    int nb_wf1 = 0;

    int nb_wf0_max = 10;
    int nb_wf1_max = 10;

    TFile *f = new TFile(output.c_str(), "RECREATE");
    TDirectory * pulse_dir = f->mkdir("pulse");
    TDirectory * wf0_dir = pulse_dir->mkdir("wfType_0");
    TDirectory * wf0_outlier_dir = wf0_dir->mkdir("outlier");
    TDirectory * wf1_dir = pulse_dir->mkdir("wfType_1");
    TDirectory * wf0_truncated_dir = pulse_dir->mkdir("wfType_0_truncated");
    TDirectory * wf0_functional_dir = pulse_dir->mkdir("wfType_0_functional");
    

    /// --- Loop over files
    int nfiles = 0;
    for (std::string filename : filenames) {
        /// --- Open files
        nfiles++;
        printf("> Open filename %d/%d: %s\n", nfiles, (int) filenames.size(), filename.c_str());
        hipo::reader  reader(filename.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);
        
        /// --- Bank definition
        hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
        hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
        hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
        hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
        hipo::bank  recBank(factory.getSchema("REC::Particle"));
        hipo::event event;
        
        /// --- Event counters
        long unsigned int nevents =0;

        bool canStop = false;

        /// --- Loop over events
        while( reader.next()){

            //if (canStop) break;
            
            /// --- display progress Bar
            nevents++;
            if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
                futils::progressBar(100.0*nevents/reader.getEntries());
            }

            /// --- load bank content for this event
            reader.read(event);
            event.getStructure(adcBank);
            event.getStructure(wfBank);
            event.getStructure(trackBank);
            event.getStructure(hitBank);
            event.getStructure(recBank);

            //bool FT_flag = false;
            
                            for (int j = 0; j < hitBank.getRows(); j++) {

                                    residual = hitBank.getDouble("residual", j);
                                    time = hitBank.getDouble("time", j);
                                    
                                    
                                    int adcRow = hitBank.getShort("id", j) - 1;
                                    adc = hitBank.getInt("adc", j);
                                    raw_adc = adcBank.getInt("ADC", adcRow);
                                    tot = adcBank.getFloat("timeOverThreshold", adcRow);
                                    leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", adcRow);
                                    wfType = adcBank.getShort("wfType", adcRow);

                                    //printf("rawdADC %d , calibrated ADC : %d \n", rawAdc, adc);

                                    Histos->H1_hit_residual->Fill(residual);
                                    Histos->H1_hit_adc->Fill(raw_adc);
                                    Histos->H1_hit_time->Fill(time);
                                    Histos->H1_hit_leadingEdgeTime->Fill(leadingEdgeTime);
                                    Histos->H1_hit_tot->Fill(tot);

                                    Histos->H2_hit_adc_vs_time->Fill(raw_adc, time);
                                    Histos->H2_hit_adc_vs_leadingEdgeTime->Fill(raw_adc, leadingEdgeTime);
                                    Histos->H2_hit_adc_vs_tot->Fill(raw_adc, tot);
                                    Histos->H2_hit_time_vs_tot->Fill(time, tot);

                                    ///////////////////////////
                                    // pulse analysis
                                    ///////////////////////////
                                    TGraph * gr = new TGraph();
                                    for (int bin  = 0; bin < 20; bin++) {
                                        double s = wfBank.getShort(TString::Format("s%d",bin+1).Data(), adcRow);
                                        gr->AddPoint(bin*48.0, s);
                                    }
                                    TGraph* gr_slope = getSlopes(gr);
                                    // extract max slope
                                    int N = gr->GetN();
                                    double ymax = 0;
                                    double xmax = 0;
                                    for (int i = 0; i < N; i++) {
                                        double x = gr_slope->GetPointX(i);
                                        double y = gr_slope->GetPointY(i);
                                        if (y > ymax) {
                                            ymax = y;
                                            xmax = x;
                                        }
                                    }
                                    slope_max = ymax;
                                    time_slope_max = xmax;
                                    Histos->H2_hit_slope_vs_adc_allTypes->Fill(raw_adc, slope_max);
                                    
                                    // Fill tree
                                    tree->Fill();

                                    // saturated signal
                                    if (wfType == 1) {

                                        //if (raw_adc > 3600)
                                            Histos->H1_hit_residual_wfType1->Fill(residual);

                                        // int smax = 0;
                                        // TGraph * gr_cut = new TGraph();
                                        // // get y max
                                        // for (int i  = 0; i < N; i++) {
                                        //     double y = gr->GetPointY(i);
                                        //     if (y > smax) smax = y;
                                        // }
                                        // // fill points and exclused plateau
                                        // for (int i  = 0; i < N; i++) {
                                        //     double x = gr->GetPointX(i);
                                        //     double y = gr->GetPointY(i);
                                        //     if (y < 0.95*smax && i > 4)
                                        //         gr_cut->AddPoint(x, y);
                                        // }
                                        nb_wf1++;
                                        
                                        // // fit gr_cut
                                        // double timeMax = adcBank.getFloat("time", adcRow);
                                        // double integral = adcBank.getInt("integral", adcRow);
                                        // double ped = adcBank.getFloat("ped", adcRow);
                                        // TF1* fitFcn = new TF1("fitFcn", fitFunction, 0, 1000, 4);
                                        // fitFcn->SetParameter(0, timeMax); // mpv --> timeMax
                                        // fitFcn->SetParameter(1, tot/4.017); // sigma --> ToT
                                        // fitFcn->SetParameter(2, integral); // norm --> integral
                                        // fitFcn->SetParameter(3, ped); // pedestal

                                        // //fitFcn->SetParLimits(0, 0, 1000);
                                        // //fitFcn->SetParLimits(1, 0, 1.5*tot);
                                        // //fitFcn->SetParLimits(2, 0, 5*integral);
                                        // //fitFcn->SetParLimits(3, -1e10, 2000);
                                        
                                        // if (nb_wf1 <= nb_wf1_max) {
                                        //     gr_cut->Fit(fitFcn, "RW");
                                        // }
                                        // else {
                                        //     gr_cut->Fit(fitFcn, "RQW");
                                        // }
                                            
                                        // double mpv = fitFcn->GetParameter(0);
                                        // //double sigma = fitFcn->GetParameter(1);
                                        // //double norm = fitFcn->GetParameter(2);
                                        // double constant = fitFcn->GetParameter(3);
                                        // double new_adc = fitFcn->GetMaximum() - constant;
                                        // //printf("mpv : %lf \n", mpv);
                                        // Int_t oldLevel = gErrorIgnoreLevel; // ignore non fatal error
                                        // gErrorIgnoreLevel = kFatal;
                                        // double t1 = fitFcn->GetX(constant + 0.5*new_adc, 0, mpv);
                                        // double t2 = fitFcn->GetX(constant + 0.5*new_adc, mpv, 1000);
                                        // gErrorIgnoreLevel = oldLevel; // end ignore non fatal error
                                        // double new_tot = t2-t1;

                                        // save some events
                                        if (nb_wf1 < nb_wf1_max) { // pulse
                                            wf1_dir->cd();
                                            { // original pulse
                                                TCanvas* c = new TCanvas();
                                                //gr->SetTitle(TString::Format("ADC = %d, ToT = %.2lf --> fitted ADC = %.2lf, ToT = %.2lf, chi2 = %lf; time (ns); ADC", raw_adc, tot, new_adc, new_tot, fitFcn->GetChisquare()).Data());
                                                gr->SetLineStyle(10);
                                                gr->Draw("APL");
                                                // gr_cut->SetMarkerColor(kRed);
                                                // gr_cut->SetMarkerStyle(20);
                                                // gr_cut->Draw("PL");
                                                c->Write(TString::Format("pulse_%d", nb_wf1).Data());
                                            }
                                            { // pulse slope
                                                wf1_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                TGraph* gr_slope = getSlopes(gr);
                                                gr_slope->Draw("APL");
                                                c->Write(TString::Format("pulse_%d_slope", nb_wf1).Data());
                                            }
                                                
                                        }

                                    } // end wfType 1

                                    else if (wfType == 0) {

                                        nb_wf0++;
                                        double predicted_adc = -28.753520 + 191.791679*slope_max;
                                        double delta_adc = raw_adc - predicted_adc;

                                        Histos->H2_hit_slope_vs_adc_wfType0->Fill(raw_adc, slope_max);
                                        Histos->H1_diff_adc_wfType0->Fill(delta_adc);

                                        //printf("raw_adc : %d , predicted_adc : %lf , slope_max : %lf \n", raw_adc, predicted_adc, slope_max);

                                        if (nb_wf0 < nb_wf0_max) {
                                            { // original pulse
                                                wf0_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                gr->SetMarkerColor(kRed);
                                                gr->SetMarkerStyle(20);
                                                gr->SetTitle(TString::Format("ADC = %d , time = %.2lf , ToT = %.2lf; time (ns); ADC", raw_adc, time, tot).Data());
                                                gr->Draw("APL");
                                                c->Write(TString::Format("pulse_%d", nb_wf0).Data());
                                            }

                                            { // pulse slope
                                                wf0_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                gr_slope->SetMarkerColor(kRed);
                                                gr_slope->SetMarkerStyle(20);
                                                gr_slope->Draw("APL");
                                                c->Write(TString::Format("pulse_%d_slope", nb_wf0).Data());
                                            }

                                        }

                                        if (nb_wf0_outlier < nb_wf0_max && fabs(delta_adc) > 200) {
                                            nb_wf0_outlier++;
                                            { // original pulse
                                                wf0_outlier_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                gr->SetMarkerColor(kRed);
                                                gr->SetMarkerStyle(20);
                                                gr->SetTitle(TString::Format("ADC = %d , time = %.2lf , ToT = %.2lf; time (ns); ADC", raw_adc, time, tot).Data());
                                                gr->Draw("APL");
                                                c->Write(TString::Format("pulse_%d", nb_wf0).Data());
                                            }

                                            { // pulse slope
                                                wf0_outlier_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                gr_slope->SetMarkerColor(kRed);
                                                gr_slope->SetMarkerStyle(20);
                                                gr_slope->SetTitle(TString::Format("Derivate of the pulse : raw ADC = %d , predicted ADC = %.2lf ", raw_adc, predicted_adc).Data());
                                                gr_slope->Draw("APL");
                                                c->Write(TString::Format("pulse_%d_slope", nb_wf0).Data());
                                            }

                                        }

                                        if (raw_adc > 1800 && raw_adc < 3500) {
                                            int N = gr->GetN();
                                            //double ADC_LIMIT = 2000;
                                            //double half_amp = ADC_LIMIT/2;
                                            double threshold = 1800;
                                            // simulate a saturated signal
                                            // TGraph* gr_truncated = new TGraph();
                                            // for (int i = 0; i < N; i++) {
                                            //     double x = gr->GetPointX(i);
                                            //     double y = gr->GetPointY(i);
                                            //     gr_truncated->AddPoint(x, std::min(y, ADC_LIMIT));
                                            // }
                                            // extract time at half amplitude
                                            double t1 = -99;
                                            double t2 = -99;
                                            double baseline = (gr->GetPointY(0) + gr->GetPointY(1) + gr->GetPointY(2) + gr->GetPointY(3))/4;
                                            for (int i = 0; i < N-1; i++) {
                                                double x1 = gr->GetPointX(i);
                                                double y1 = gr->GetPointY(i) - baseline;
                                                double x2 = gr->GetPointX(i+1);
                                                double y2 = gr->GetPointY(i+1) - baseline;
                                                double a = (y1-y2)/(x1-x2);
                                                if (y1 < threshold && y2 > threshold) {
                                                    t1 = x1 + (threshold-y1)/a;
                                                }
                                                if (y1 > threshold && y2 < threshold && y2 > 1) { // make sure the second point is not zero
                                                    t2 = x1 + (threshold-y1)/a;
                                                }
                                            }
                                            if (t1 < 0 || t2 < 0) continue;
                                            double new_tot = t2-t1;
                                            Histos->H2_hit_new_tot_vs_adc->Fill(raw_adc, new_tot);
                                            Histos->H2_hit_slope_vs_adc_comp_new_tot->Fill(raw_adc, slope_max);
                                            if (nb_wf0_truncated < nb_wf0_max) {
                                                nb_wf0_truncated++;
                                                wf0_truncated_dir->cd();
                                                TCanvas* c = new TCanvas();
                                                gr->SetMarkerColor(kRed);
                                                gr->SetMarkerStyle(20);
                                                gr->SetTitle(TString::Format("ADC = %d , time = %.2lf , ToT = %.2lf; time (ns); ADC", raw_adc, time, tot).Data());
                                                gr->Draw("APL");
                                                // gr_truncated->SetLineColor(kBlue);
                                                // gr_truncated->SetLineStyle(10);
                                                // gr_truncated->Draw("L same");
                                                c->Write(TString::Format("pulse_%d", nb_wf0).Data());
                                            }
                                        }

                                        // study the functional
                                        if (raw_adc > 2000 && raw_adc < 2200) {
                                            nb_wf0_functional++;
                                            if (nb_wf0_functional < 30) {
                                                { // pulse
                                                    wf0_functional_dir->cd();
                                                    TCanvas* c = new TCanvas();
                                                    gr->SetMarkerColor(kRed);
                                                    gr->SetMarkerStyle(20);
                                                    gr->SetTitle(TString::Format("ADC = %d , time = %.2lf , ToT = %.2lf; time (ns); ADC", raw_adc, time, tot).Data());
                                                    gr->Draw("APL");
                                                    c->Write(TString::Format("pulse_%d", nb_wf0_functional).Data());
                                                    
                                                    //////////////
                                                    /// Fit
                                                    //////////////
                                                    int Neff = gr->GetN();
                                                    for (int i = gr->GetN()-1; i >= 0; i--) {
                                                        if (gr->GetPointY(i) > 1) { // I want bigger than 0, but 1 if enough
                                                            Neff = i+1;
                                                            break;
                                                        }
                                                    }
                                                    double baseline = (gr->GetPointY(0) + gr->GetPointY(1) + gr->GetPointY(2) + gr->GetPointY(3))/4;
                                                    // create histogram
                                                    TH1D* h = new TH1D(TString::Format("hist_pulse_%d", nb_wf0_functional).Data(), "Landau conv. Gauss", Neff, -24, 48*Neff-24);
                                                    for (int i = 0; i < Neff; i++) {
                                                        h->Fill(gr->GetPointX(i), std::max(gr->GetPointY(i)-baseline,0.));
                                                    }

                                                    double peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
                                                    double sigma = h->GetStdDev();
                                                    double integral = h->Integral();
                                                    

                                                       // Setting fit range and start values
                                                    double fr[2];
                                                    double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
                                                    // fr[0]=0.3*h->GetMean();
                                                    // fr[1]=3.0*h->GetMean();
                                                    fr[0]=peak-sigma;
                                                    fr[1]=peak+1.5*sigma;
                                                    // fr[0]=peak-1.5*sigma;
                                                    // fr[1]=peak+0.05*sigma;
                                                    // fr[0]=peak-2*sigma;;
                                                    // fr[1]=0.7*peak;
                                                    
                                                    pllo[0]=0.6*tot/4.017; pllo[1]=peak-2*sigma; pllo[2]=1.0; pllo[3]=0.4;
                                                    plhi[0]=2*tot/4.017; plhi[1]=peak+2*sigma; plhi[2]=2000000.0; plhi[3]=15.0;
                                                    sv[0]=tot/4.017; sv[1]=peak; sv[2]=integral; sv[3]=10;

                                                    double chisqr;
                                                    int    ndf;
                                                    TF1 *fitsnr = langaufit(h,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

                                                    //h->Write(TString::Format("hist_pulse_%d", nb_wf0_functional).Data());
                                                    TCanvas* c2 = new TCanvas();
                                                    h->Draw();
                                                    TF1* fitsnr_copy = (TF1*) fitsnr->Clone(TString::Format("%s_copy", fitsnr->GetName()).Data());
                                                    fitsnr->SetRange(0,1000);
                                                    fitsnr->SetLineStyle(2);
                                                    fitsnr->Draw("l same");
                                                    fitsnr_copy->Draw("l same");
                                                    c2->Write(TString::Format("hist_pulse_%d", nb_wf0_functional).Data());

                                                }
                                                { // pulse slope
                                                    wf0_functional_dir->cd();
                                                    TCanvas* c = new TCanvas();
                                                    TGraph* gr_slope = getSlopes(gr);
                                                    gr_slope->Draw("APL");
                                                    c->Write(TString::Format("pulse_%d_slope", nb_wf0_functional).Data());
                                                }
                                            } 
                                            // else {
                                            //     //canStop = true;
                                            // }
                                        }

                                    } // end wfType 0

                                } // loop over hits

        } // end loop over events for this file


    } // loop over files

    /// --- Output
    
    fitGaus(Histos->H1_hit_residual_wfType1);

    Histos->SaveIn(f);


    /// --- post processing

    { // adc versus slope
        if (f->cd("time_versus_adc")) {
            f->mkdir("time_versus_adc/fits_slope_vs_adc");
        }
        TH2D* h2 = (TH2D*) Histos->H2_hit_slope_vs_adc_wfType0->Clone("adc_vs_slopemax");
        TGraphErrors* gr = new TGraphErrors();
        int nbins = h2->GetXaxis()->GetNbins();
        int step = 5;
        int bin = 4;
        while (bin < nbins && bin < nbins-3*step) {
            int binSup = std::min(bin+step-1, nbins);
            TH1D* h_tmp = h2->ProjectionY(TString::Format("h_tmp_%d-%d", bin, binSup).Data(), bin, binSup);
            // fit
            h_tmp->Fit("gaus", "RQ", "", h_tmp->GetMean() - 2*h_tmp->GetStdDev(), h_tmp->GetMean() + 2*h_tmp->GetStdDev());
            TF1 *gaus = h_tmp->GetFunction("gaus");
            double mean = gaus->GetParameter("Mean");
            double sigma = gaus->GetParameter("Sigma");
            // data
            double xinf = h2->GetXaxis()->GetBinCenter(bin);
            double xsup = h2->GetXaxis()->GetBinCenter(binSup);
            double xval = 0.5*(xinf+xsup);
            //double yval = h_tmp->GetMaximum();
            // save
            gr->AddPointError(xval, mean, 0, sigma);
            if (f->cd("time_versus_adc/fits_slope_vs_adc")) 
                h_tmp->Write(h_tmp->GetName());
            // go to next bin
            bin += step;
        }

        if (f->cd("time_versus_adc")) {
            TCanvas* c = new TCanvas();
            h2->Draw("colz");
            h2->SetStats(0);
            gr->SetMarkerColor(kBlack);
            gr->SetMarkerStyle(21);
            gr->SetMarkerSize(2);
            gr->SetLineWidth(3);
            gr->SetLineColor(kBlack);
            gr->Draw("same pe");
            // fit
            gr->Fit("pol1","R","", 0, 4000);
            TF1* fn = gr->GetFunction("pol1");
            fn->SetLineWidth(3);
            TText text;
            text.SetTextSize(0.035);
            double p0 = fn->GetParameter(0);
            double p1 = fn->GetParameter(1);
            text.DrawText(250, 18, TString::Format("slope = %lf + %lf*ADC", p0, p1).Data());
            text.DrawText(250, 16, TString::Format("ADC = %lf + %lf*slope", -p0/p1, 1/p1).Data());
            c->Write("hit_slope_vs_adc_wfType0_fitted");
        }
    }

    TGraph * gr_res_tot = new TGraph();
    TGraph * gr_res_slope = new TGraph();

    {
        // delta ADC versus tot for 1800 < ADC < 3500
        if (f->cd("time_versus_adc")) {
            f->mkdir("time_versus_adc/fits_delta_adc_vs_tot");
        }
        TH2D* h2 = (TH2D*) Histos->H2_hit_new_tot_vs_adc->Clone("new_tot_delta_adc");
        TGraphErrors* gr = new TGraphErrors();
        int nbins = h2->GetYaxis()->GetNbins();
        int step = 3;
        int bin = 8;
        while (bin < 32) {
            int binSup = std::min(bin+step-1, nbins);
            double yinf = h2->GetYaxis()->GetBinCenter(bin);
            double ysup = h2->GetYaxis()->GetBinCenter(binSup);
            double yval = 0.5*(yinf+ysup);
            TH1D* h_tmp = h2->ProjectionX(TString::Format("h_tmp_%d-%d", bin, binSup).Data(), bin, binSup);
            // fit
            double peak = h_tmp->GetXaxis()->GetBinCenter(h_tmp->GetMaximumBin());
            double dev = h_tmp->GetStdDev();
            h_tmp->Fit("gaus", "RQ", "", peak-dev, peak+dev);
            TF1 *gaus = h_tmp->GetFunction("gaus");
            double mean = gaus->GetParameter("Mean");
            double sigma = gaus->GetParameter("Sigma");
            // save
            //gr->AddPointError(mean, yval, sigma, 0);
            gr->AddPointError(peak, yval, dev, 0);
            if (f->cd("time_versus_adc/fits_delta_adc_vs_tot")) 
                h_tmp->Write(h_tmp->GetName());
            // go to next bin
            bin += step;
        }
        gr->SetTitle("ADC resolution versus ToT; ToT (ns); ADC");
        gr->Write("gr_delta_adc_vs_tot");
        int N = gr->GetN();
        TGraph* grerr = new TGraph();
        for (int i = 0; i < N; i++) {
            double x = gr->GetPointY(i);
            double y = gr->GetErrorX(i);
            grerr->AddPoint(x,y);
            gr_res_tot->AddPoint(gr->GetPointX(i), gr->GetErrorX(i));
        }
        grerr->SetTitle("ADC resolution; ToT; ADC resolution");
        grerr->Write("error_graph_delta_adc_vs_tot");

        // canvas
        TCanvas* c = new TCanvas();
        h2->Draw("colz");
        h2->SetStats(0);
        gr->SetMarkerColor(kBlack);
        gr->SetMarkerStyle(21);
        gr->SetMarkerSize(2);
        gr->SetLineWidth(3);
        gr->SetLineColor(kBlack);
        gr->Draw("same pe");
        c->Write("hit_new_tot_vs_adc_fitted");
    }

    {
        // delta ADC versus tot for 1800 < ADC < 3500
        if (f->cd("time_versus_adc")) {
            f->mkdir("time_versus_adc/fits_delta_adc_vs_slope");
        }
        TH2D* h2 = (TH2D*) Histos->H2_hit_slope_vs_adc_comp_new_tot->Clone("slope_delta_adc_comp_new_tot");
        TGraphErrors* gr = new TGraphErrors();
        int nbins = h2->GetYaxis()->GetNbins();
        int step = 3;
        int bin = 17;
        while (bin < 46) {
            int binSup = std::min(bin+step-1, nbins);
            double yinf = h2->GetYaxis()->GetBinCenter(bin);
            double ysup = h2->GetYaxis()->GetBinCenter(binSup);
            double yval = 0.5*(yinf+ysup);
            TH1D* h_tmp = h2->ProjectionX(TString::Format("h_tmp_%d-%d", bin, binSup).Data(), bin, binSup);
            // fit
            double peak = h_tmp->GetXaxis()->GetBinCenter(h_tmp->GetMaximumBin());
            double dev = h_tmp->GetStdDev();
            h_tmp->Fit("gaus", "RQ", "", peak-dev, peak+dev);
            TF1 *gaus = h_tmp->GetFunction("gaus");
            double mean = gaus->GetParameter("Mean");
            double sigma = gaus->GetParameter("Sigma");

            // save
            //gr->AddPointError(mean, yval, sigma, 0);
            gr->AddPointError(peak, yval, dev, 0);
            if (f->cd("time_versus_adc/fits_delta_adc_vs_slope")) 
                h_tmp->Write(h_tmp->GetName());
            // go to next bin
            bin += step;
        }
        gr->SetTitle("ADC resolution versus Slope; ADC; Slope (ADC/ns)");
        gr->Write("gr_delta_adc_vs_slope");
        int N = gr->GetN();
        TGraph* grerr = new TGraph();
        for (int i = 0; i < N; i++) {
            double x = gr->GetPointY(i);
            double y = gr->GetErrorX(i);
            grerr->AddPoint(x,y);
            gr_res_slope->AddPoint(gr->GetPointX(i), gr->GetErrorX(i));
        }
        grerr->SetTitle("ADC resolution; Slope; ADC resolution");
        grerr->Write("error_graph_delta_adc_vs_slope");

        // canvas
        TCanvas* c = new TCanvas();
        h2->Draw("colz");
        h2->SetStats(0);
        gr->SetMarkerColor(kBlack);
        gr->SetMarkerStyle(21);
        gr->SetMarkerSize(2);
        gr->SetLineWidth(3);
        gr->SetLineColor(kBlack);
        gr->Draw("same pe");
        c->Write("hit_slope_vs_adc_comp_new_tot");
    }

    // compare error resolution for tot and slope
    {
        TCanvas* c = new TCanvas();
        TMultiGraph* mg = new TMultiGraph();
        TLegend* legend = new TLegend();

        gr_res_slope->SetLineColor(kRed);
        gr_res_slope->SetLineWidth(2);
        //gr_res_slope->Draw("APL");
        mg->Add(gr_res_slope, "pl");

        gr_res_tot->SetLineColor(kBlue);
        gr_res_tot->SetLineWidth(2);
        //gr_res_tot->Draw("PL same");
        mg->Add(gr_res_tot, "pl");

        mg->SetTitle("ADC resolution for each method; ADC; ADC resolution");
        mg->Draw("APL");

        legend->AddEntry(gr_res_slope, "slope", "lpf");
        legend->AddEntry(gr_res_tot, "ToT", "lpf");
        legend->Draw();

        if (f->cd("time_versus_adc")) {
            c->Write("comparison_slope_vs_tot_ADC_resolution");
        }
    }




    

    

    TDirectory * tree_dir = f->mkdir("data");
    tree_dir->cd();
    tree->Write();


    f->Close();
    printf("> Output file created : %s\n", output.c_str());
    
    // end of the program
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("\033[1m * time elapsed : %lf seconds\033[0m\n", elapsed.count());
    return 0;
}


// ./lookat10gev.exe -i /volatile/clas12/rg-l/production/p0v10_3_alert/mon/recon/021459/rec_rec_clas_021459.evio.00000.hipo -o ./output/adc-versus-time-p0v10_3.root