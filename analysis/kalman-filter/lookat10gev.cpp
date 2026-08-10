/***********************************************
 * Look at 10.6 GeV data
 *
 * @author Felix Touchte Codjo
 * @date April 02, 2026
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

#include "lookat10gev.h"
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

    tree->Branch("residual", &residual, "residual/D");
    tree->Branch("time", &time, "time/D");
    tree->Branch("leadingEdgeTime", &leadingEdgeTime, "leadingEdgeTime/D");
    tree->Branch("tot", &tot, "tot/D");
    tree->Branch("adc", &adc, "adc/I");
    tree->Branch("raw_adc", &raw_adc, "raw_adc/I");
    tree->Branch("wfType", &wfType, "wfType/I");
    tree->Branch("slope_max", &slope_max, "slope_max/D");


    int nb_wf0 = 0;
    int nb_wf1 = 0;

    int nb_wf0_max = 10;
    int nb_wf1_max = 10;

    TFile *f = new TFile(output.c_str(), "RECREATE");
    TDirectory * pulse_dir = f->mkdir("pulse");
    TDirectory * wf0_dir = pulse_dir->mkdir("wfType_0");
    TDirectory * wf1_dir = pulse_dir->mkdir("wfType_1");
    

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

        /// --- Loop over events
        while( reader.next()){
            
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
            
            /// --- Loop over electron
            for (int i = 0; i < recBank.getRows(); i++) {
                // Look for electrons
                if (recBank.getInt("pid", i) == 11) {
                    //int status = trackBank.getShort("status", i);
                    double vze = recBank.getFloat("vz", i) * Units::cm;
                    //int status = recBank.getShort("status", i);
                    //printf("%f  ", recBank.getFloat("vz", i));
                    // if (fabs(vze+3) < 1e-6) {
                    //     FT_flag = true;
                    //     recBank.show();
                    // }
                    Histos->H1_electron_vz_nocut->Fill(vze);
                    // if (abs(status)/1000 == 1) { // FT electron
                    //     FT_flag = true;
                    //     //recBank.show();
                    // }
                    // if (FT_flag) continue;

                    {
                        double px = recBank.getFloat("px", i) * Units::GeV;
                        double py = recBank.getFloat("py", i) * Units::GeV;
                        double pz = recBank.getFloat("pz", i) * Units::GeV;
                        double p = sqrt(px*px + py*py + pz*pz);
                        double theta = acos(pz/p);
                        double phi = atan2(py, px);

                        Histos->H1_electron_p->Fill(p);
                        Histos->H1_electron_theta->Fill(theta / Units::deg); // read per degree to convert in degree
                        Histos->H1_electron_phi->Fill(phi > 0 ? phi / Units::deg : 360 + phi / Units::deg);
                        Histos->H1_electron_vz->Fill(vze);
                    }
                    
                    /// --- Loop over track
                    for (int i = 0; i < trackBank.getRows(); i++) {
                        
                        int nhits = trackBank.getInt("n_hits", i);
                        double chi2 = trackBank.getFloat("chi2", i);

                        if (nhits >= 7 && chi2 < 8) {
                        
                            double vz = trackBank.getFloat("z", i) * Units::mm;
                            double px = trackBank.getFloat("px", i) * Units::MeV;
                            double py = trackBank.getFloat("py", i) * Units::MeV;
                            double pz = trackBank.getFloat("pz", i) * Units::MeV;
                            double dEdx = trackBank.getFloat("dEdx", i) * Units::MeV/Units::mm;
                            

                            

                            double p = sqrt(px*px + py*py + pz*pz);
                            double theta = acos(pz/p);
                            double phi = atan2(py, px);

                            Histos->H1_track_p->Fill(p);
                            Histos->H1_track_theta->Fill(theta / Units::deg); // read per degree to convert in degree
                            Histos->H1_track_phi->Fill(phi > 0 ? phi / Units::deg : 360 + phi / Units::deg);
                            Histos->H1_track_nhits->Fill(nhits);
                            Histos->H1_track_chi2->Fill(chi2);
                            Histos->H1_track_vz->Fill(vz);
                            Histos->H2_track_corr_p_dEdx->Fill(p, dEdx / (Units::MeV/Units::mm));

                            Histos->H1_delta_vz->Fill(vze - vz);

                            // Loop over hit for this track
                            int trackid = trackBank.getInt("trackid", i);
                            for (int j = 0; j < hitBank.getRows(); j++) {
                                if (hitBank.getInt("trackid", j) == trackid) {
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

                                    Histos->H2_hit_adc_vs_time->Fill(adc, time);
                                    Histos->H2_hit_adc_vs_leadingEdgeTime->Fill(adc, leadingEdgeTime);
                                    Histos->H2_hit_adc_vs_tot->Fill(adc, tot);
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
                                    //double xmax = 0;
                                    for (int i = 0; i < N; i++) {
                                        //double x = gr_slope->GetPointX(i);
                                        double y = gr_slope->GetPointY(i);
                                        if (y > ymax) {
                                            ymax = y;
                                            //xmax = x;
                                        }
                                    }
                                    slope_max = ymax;
                                    Histos->H2_hit_slope_vs_adc_allTypes->Fill(raw_adc, slope_max);
                                    
                                    // Fill tree
                                    tree->Fill();

                                    // saturated signal
                                    if (wfType == 1) {
                                        int smax = 0;
                                        TGraph * gr_cut = new TGraph();
                                        // get y max
                                        for (int i  = 0; i < N; i++) {
                                            double y = gr->GetPointY(i);
                                            if (y > smax) smax = y;
                                        }
                                        // fill points and exclused plateau
                                        for (int i  = 0; i < N; i++) {
                                            double x = gr->GetPointX(i);
                                            double y = gr->GetPointY(i);
                                            if (y < 0.95*smax && i > 4)
                                                gr_cut->AddPoint(x, y);
                                        }
                                        nb_wf1++;
                                        
                                        // fit gr_cut
                                        double timeMax = adcBank.getFloat("time", adcRow);
                                        double integral = adcBank.getInt("integral", adcRow);
                                        double ped = adcBank.getFloat("ped", adcRow);
                                        TF1* fitFcn = new TF1("fitFcn", fitFunction, 0, 1000, 4);
                                        fitFcn->SetParameter(0, timeMax); // mpv --> timeMax
                                        fitFcn->SetParameter(1, tot/4.017); // sigma --> ToT
                                        fitFcn->SetParameter(2, integral); // norm --> integral
                                        fitFcn->SetParameter(3, ped); // pedestal

                                        //fitFcn->SetParLimits(0, 0, 1000);
                                        //fitFcn->SetParLimits(1, 0, 1.5*tot);
                                        //fitFcn->SetParLimits(2, 0, 5*integral);
                                        //fitFcn->SetParLimits(3, -1e10, 2000);
                                        
                                        if (nb_wf1 <= nb_wf1_max) {
                                            gr_cut->Fit(fitFcn, "RW");
                                        }
                                        else {
                                            gr_cut->Fit(fitFcn, "RQW");
                                        }
                                            
                                        double mpv = fitFcn->GetParameter(0);
                                        //double sigma = fitFcn->GetParameter(1);
                                        //double norm = fitFcn->GetParameter(2);
                                        double constant = fitFcn->GetParameter(3);
                                        double new_adc = fitFcn->GetMaximum() - constant;
                                        //printf("mpv : %lf \n", mpv);
                                        Int_t oldLevel = gErrorIgnoreLevel; // ignore non fatal error
                                        gErrorIgnoreLevel = kFatal;
                                        double t1 = fitFcn->GetX(constant + 0.5*new_adc, 0, mpv);
                                        double t2 = fitFcn->GetX(constant + 0.5*new_adc, mpv, 1000);
                                        gErrorIgnoreLevel = oldLevel; // end ignore non fatal error
                                        double new_tot = t2-t1;

                                        // save some events
                                        if (nb_wf1 < nb_wf1_max) { // pulse
                                            wf1_dir->cd();
                                            { // original pulse
                                                TCanvas* c = new TCanvas();
                                                gr->SetTitle(TString::Format("ADC = %d, ToT = %.2lf --> fitted ADC = %.2lf, ToT = %.2lf, chi2 = %lf; time (ns); ADC", raw_adc, tot, new_adc, new_tot, fitFcn->GetChisquare()).Data());
                                                gr->SetLineStyle(10);
                                                gr->Draw("AL");
                                                gr_cut->SetMarkerColor(kRed);
                                                gr_cut->SetMarkerStyle(20);
                                                gr_cut->Draw("PL");
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
                                        Histos->H2_hit_slope_vs_adc_wfType0->Fill(raw_adc, slope_max);

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
                                        
                                    } // end wfType 0

                                    

                                }
                            } // loop over track hits

                        } // track selection

                    } // loop over tracks 
                    
                } // electron selection
            }

        } // end loop over events for this file


    } // loop over files

    /// --- Output
    
    Histos->SaveIn(f);


    /// --- post processing
    if (f->cd("time_versus_adc")) {
        f->mkdir("time_versus_adc/fits");
    }
    TH2D* h2 = (TH2D*) Histos->H2_hit_slope_vs_adc_wfType0->Clone("adc_vs_slopemax");
    TGraphErrors* gr = new TGraphErrors();
    int nbins = h2->GetXaxis()->GetNbins();
    int step = 5;
    int bin = 4;
    while (bin < nbins && bin < nbins-3*step) {
        int binSup = std::min(bin+step, nbins);
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
        if (f->cd("time_versus_adc/fits")) 
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
        text.DrawText(250, 16, TString::Format("ADC = %lf + %lf*slope", -p1/p0, 1/p0).Data());
        c->Write("hit_slope_vs_adc_wfType0_fitted");
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