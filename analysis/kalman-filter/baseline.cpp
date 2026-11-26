/***********************************************
 * Study the baseline of the waveforms
 *
 * @author Felix Touchte Codjo
 * @date November 26, 2026
 * 
 * ********************************************/

#include <cstdio>
#include <cmath>
#include <random>

#include <vector>
#include <string>
#include <chrono>

#include "fOptions.h"
#include "Units.h"

#include "../hipo4/reader.h"
#include "../hipo4/writer.h"

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
#include "TLine.h"


/**
 * @brief Return an angle between 0 and 2*pi
 * 
 * @param x 
 * @return double 
 */
double mod2Pi(double x) {
    double r = atan2(sin(x), cos(x)); // is an angle between -pi and pi
    if (r >= 0) {
        return r;
    } else {
        return r + 2*M_PI;
    }
}

/**
 * @brief Return an angle between -pi and +pi
 * 
 * equals atan2(sin(x), cos(x))
 * 
 * @param x 
 * @return double 
 */
double modPi(double x) {
    return atan2(sin(x), cos(x));
}

/**
 * @brief Progress bar
 * 
 * @param state is a number below bar_length
 * @param bar_length is 100 per default
 */
void progressBar(int state, int bar_length = 100) { // state is a number between 0 and 100
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

std::vector<double> linear_regression(const std::vector<double> & x, const std::vector<double> & y) {
    if (x.size() != y.size()) return {0,0};
    int N = x.size();

    double sum_x = 0;
    double sum_y = 0;
    double sum_x2 = 0;
    double sum_y2 = 0;
    double sum_xy = 0;

    for (int i = 0; i < N; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i]*x[i];
        sum_y2 += y[i]*y[i];
        sum_xy += x[i]*y[i];
    }

    std::vector<double> b(2);
    double denum = N*sum_x2 - sum_x*sum_x;
    b[0] = (sum_y*sum_x2 - sum_x*sum_xy)/denum;
    b[1] = (N*sum_xy - sum_x*sum_y)/denum;

    return b;
}

void test_linear_regression() {
    int N = 100;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution d{0.0, 1.0};
    std::vector<double> x(N), y(N);
    for (int i = 0; i < N; i++) {
        x[i] = sin(i);
        y[i] = 5 + d(gen) + 4*x[i];
    }
    std::vector<double> b = linear_regression(x,y);
    printf("b[0] = %lf\n", b[0]);
    printf("b[1] = %lf\n", b[1]);

    TCanvas* c = new TCanvas();
    TGraph* gr = new TGraph();
    for (int i = 0; i < N; i++) {
        gr->AddPoint(x[i], y[i]);
    }
    TF1* f1 = new TF1("fit", "[0] + [1]*x", -1, 1);
    f1->SetParameter(0, b[0]);
    f1->SetParameter(1, b[1]);
    gr->SetMarkerStyle(8);
    gr->Draw("ap");
    f1->Draw("l same");
    c->Print("./output/test_linear_regression.pdf");
}

void process_waveform(const std::vector<double> & s, std::string name, TDirectory* dir) {
    int N = s.size();
    std::vector<double> t(2); // trailing edge
    t[0] = 0;
    t[1] = 0;
    for (int i = 0; i < N-1; i++) {
        double tmp = s[i+1]-s[i];
        if (tmp > t[1]) { // catch the maximum slope
            t[1] = tmp;
            t[0] = s[i]-t[1]*i; // y = a*(x-x0) + y0 ==> y = a*x + y0 - a*x0
        }
    }

    std::vector<double> x(4), y(4);
    for (int i = 0; i < 4; i++) {
        x[i] = i;
        y[i] = s[i];
    }
    std::vector<double> b = linear_regression(x,y); // baseline

    // Draw waveform
    dir->cd();
    TCanvas* c = new TCanvas();
    TGraph* gr = new TGraph();
    for (int i = 0; i < N; i++) {
        gr->AddPoint(i, s[i]);
    }
    gr->SetMarkerStyle(8);
    gr->SetLineColor(kBlack);
    gr->Draw("apl");

    // Draw trailing edge
    TF1* f1t = new TF1("fit", "[0] + [1]*x", 0, N);
    f1t->SetParameter(0, t[0]);
    f1t->SetParameter(1, t[1]);
    f1t->Draw("l same");

    // Draw baseline
    TF1* f1b = new TF1("fit", "[0] + [1]*x", 0, N);
    f1b->SetParameter(0, b[0]);
    f1b->SetParameter(1, b[1]);
    f1b->Draw("l same");

    c->Write(name.c_str());
}


int main(int argc, char const *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();

    //test_linear_regression();
    
    /// --- Load files from options
    fOptions OPT({"-i", "-o", "-filtering"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    TFile *f = new TFile(output.c_str(), "RECREATE");
    TDirectory* negative_slope_dir = f->mkdir("negative_slope_dir");
    TDirectory* positive_slope_dir = f->mkdir("positive_slope_dir");
    TDirectory* normal_slope_dir = f->mkdir("normal_slope_dir");


    /// --- Histograms
    TH1D* h1_leadingEdgeTime = new TH1D("leadingEdgeTime", "leadingEdgeTime; leadingEdgeTime (ns); count", 100, 0, 1000);
    TH1D* h1_baseline = new TH1D("baseline", "baseline; baseline (ADC); count", 100, 0, 700);
    TH1I* h1_wfType = new TH1I("wfType", "wfType; wfType; count", 7, 0, 7);

    TH1D* h1_baseline_constant = new TH1D("baseline_constant", "baseline constant; baseline constant (ADC); count", 100, 0, 700);
    TH1D* h1_baseline_slope = new TH1D("baseline_slope", "baseline slope; baseline slope (ADC); count", 100, -150, 150);
    TH2D* h2_baseline_corr = new TH2D("baseline_corr", "baseline slope vs constant; constant (ADC); slope (ADC)", 50, 0, 700, 50, -150, 150);
    TH2D* h2_baseline_vs_time = new TH2D("baseline_vs_time", "baseline slope vs leadingEdgeTime; leadingEdgeTime (ns); baseline slope (ADC)", 50, 0, 1000, 50, -150, 150);

    bool filtering = OPT.GetValue("-filtering") == "true";

    /// -- HIPO writer : filter events
    hipo::writer writer;

    /// --- Loop over files
    int nfiles = 0;
    for (std::string filename : filenames) {
        /// --- Open files
        nfiles++;
        printf("> Open filename %d/%d: %s\n", nfiles, (int) filenames.size(), filename.c_str());
        hipo::reader  reader(filename.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        if (nfiles == 1 && filtering) {
            writer.addDictionary(factory);
            std::string hipo_output = output;
            if (hipo_output.substr(hipo_output.size()-4, hipo_output.size()) == ".root") {
                hipo_output = hipo_output.substr(0, hipo_output.size()-4);
            }
            hipo_output += ".hipo";
            writer.open(hipo_output.c_str());
        }
        
        /// --- Bank definition
        hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
        hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
        // hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
        // hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
        // hipo::bank  recBank(factory.getSchema("REC::Particle"));
        hipo::event event;
        
        // Event counters
        long unsigned int nevents =0;

        int nb_negative_slope = 0;
        int nb_positive_slope = 0;
        int nb_normal_slope = 0;

        /// --- Loop over events
        while( reader.next()){
            
            /// --- display progress Bar
            nevents++;
            if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
                progressBar(100.0*nevents/reader.getEntries());
            }

            /// --- load bank content for this event
            reader.read(event);
            event.getStructure(adcBank);
            event.getStructure(wfBank);
            // event.getStructure(trackBank);
            // event.getStructure(hitBank);
            // event.getStructure(recBank);

            for (int row = 0; row < wfBank.getRows(); row++) {

                std::vector<double> x(4), y(4);
                for (int i = 0; i < 4; i++) {
                    char* entry = new char[5];
                    sprintf(entry, "s%d", i+1);
                    int adc = wfBank.getShort(entry, row);
                    x[i] = i;
                    y[i] = adc;
                    //printf("%s\n", entry);
                }

                std::vector<double> b = linear_regression(x,y);

                int wfType = adcBank.getShort("wfType", row);
                double leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", row);
                int ped = adcBank.getFloat("ped", row);
                //int integra

                if (wfType > 0) continue;


                h1_leadingEdgeTime->Fill(leadingEdgeTime);
                h1_wfType->Fill(wfType);
                h1_baseline->Fill(ped);
                h1_baseline_constant->Fill(b[0]);
                h1_baseline_slope->Fill(b[1]);
                h2_baseline_corr->Fill(b[0], b[1]);
                h2_baseline_vs_time->Fill(leadingEdgeTime, b[1]);

                // Fill pulse
                std::vector<double> s(20);
                for (int i = 0; i < 20; i++) {
                    char* entry = new char[9];
                    sprintf(entry, "s%d", i+1);
                    int adc = wfBank.getShort(entry, row);
                    s[i] = adc;
                }

                // Negative slope
                if (b[1] < -20) {                     
                    nb_negative_slope++;
                    if (nb_negative_slope <= 20) {
                        std::string name = "wf_";
                        name += std::to_string(nb_negative_slope);
                        process_waveform(s, name, negative_slope_dir);
                    }
                }

                // Positive slope
                else if (b[1] > 20) {                     
                    nb_positive_slope++;
                    if (nb_positive_slope <= 20) {
                        std::string name = "wf_";
                        name += std::to_string(nb_positive_slope);
                        process_waveform(s, name, positive_slope_dir);
                    }
                }

                // Negative slope
                else if (b[1] > -10 && b[1] < 10) {                     
                    nb_normal_slope++;
                    if (nb_normal_slope <= 20) {
                        std::string name = "wf_";
                        name += std::to_string(nb_normal_slope);
                        process_waveform(s, name, normal_slope_dir);
                    }
                }



            } // end loop over wf rows






        } // end loop over events

    } // end loop over files


    if (filtering) {
        writer.close();
        writer.showSummary();
    }

    // {
    //     TH1D* h1 = (TH1D*) h1_baseline_slope->Clone("h1_baseline_processed");
    //     double population = 0.6827;
    //     double alpha = 0.5*(1-population);
    //     double p[2] = {alpha, 1-alpha};
    //     double q[2];
    //     h1->GetQuantiles(2, q, p);

    //     TCanvas* c = new TCanvas();
    //     h1->Draw("hsit");
    //     double count = h1->GetMaximum();
    //     TLine* l0 = new TLine(q[0], 0, q[0], count);
    //     TLine* l1 = new TLine(q[1], 0, q[1], count);
        
    //     l0->Draw("l same");
    //     l1->Draw("l same");

    //     l0->SetLineColor(kRed);
    //     l1->SetLineColor(kRed);

    //     f->cd();
    //     c->Write("baseline_processed");
    // }
    


    // Save
    f->cd();
    h1_leadingEdgeTime->Write(h1_leadingEdgeTime->GetName());
    h1_wfType->Write(h1_wfType->GetName());
    h1_baseline->Write(h1_baseline->GetName());
    h1_baseline_constant->Write(h1_baseline_constant->GetName());
    h1_baseline_slope->Write(h1_baseline_slope->GetName());
    h2_baseline_corr->Write(h2_baseline_corr->GetName());
    h2_baseline_vs_time->Write(h2_baseline_vs_time->GetName());



    f->Close();
    printf("> Output file created : %s\n", output.c_str());
    
    // End of the program
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("\033[1m * time elapsed : %lf seconds\033[0m\n", elapsed.count());
    return 0;

}
