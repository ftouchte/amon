/********************************************************************
 * Code to check wire disposition using elastics
 * 
 * The track will go in the opposite direction of the electron. 
 * So we should a correlation between the phi of the electron and
 * phi disposition as they are numeroted in increasing phi
 * 
 * @author ftouchte
 * @date April 23, 2026
 ********************************************************************/

 
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

#include "futils.h"
#include "fOptions.h"
#include "Units.h"
#include "AhdcUtils.h"

void progressBar(int state, int bar_length = 100);
TCanvas* showCuts(TH1D* h, double lim_inf, double lim_sup);

/// --- Constants
//const double beam_energy =  10676.6 * Units::MeV;
const double beam_energy = 2.23951 * Units::GeV; // incident energy if the electron, GeV
const double electron_mass = 0.511e-3 * Units::GeV; // energy mass of electron, GeV
const double proton_mass = 938.272e-3 * Units::GeV; // energy mass of proton, GeV
const double helium_mass = 3.73 * Units::GeV; // energy mass of Helium-4, GeV
const double deuteron_mass = 1.875 * Units::GeV; // energy mass of Deuterium, GeV


int main(int argc, char const *argv[]) {

    /// --- Selection cuts
    // double chi2_min = 0.1;
    // double chi2_max = 3;
    // double vz_min = -24 * Units::cm;
    // double vz_max = 15 * Units::cm; 
    // double nhits_min = 7;
    // double delta_phi_cut = 20 * Units::deg; // deg
    double W2_min = 3.5 * Units::GeV * Units::GeV;
    double W2_max = 3.8 * Units::GeV * Units::GeV;


    /// --- start timer
    auto start = std::chrono::high_resolution_clock::now();

    /// --- Load options
    fOptions OPT({"-i", "-o", "-v", "-simu"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    /// --- Initialise histograms
    TH1D* H1_W2 = new TH1D("W2", "W^{2}; W^{2} (GeV^{2}); count", 100, 3.2, 6);

    std::vector<TH2D*> H2_corr_wire_phi;
    for (int i = 1; i <= 8; i++) {
        int layer = AhdcUtils::number2layer(i);
        int nbWires = AhdcUtils::layerNbWires(layer);
        std::string name = Form("wire_disposition_layer_%d", layer);
        std::string title = Form("Wire occupation versus Electron phi (layer %d); #phi_{e} (deg); wire", layer);
        TH2D* h2 = new TH2D(name.c_str(), title.c_str(), nbWires, 0, 360, nbWires, 1, nbWires+1);
        H2_corr_wire_phi.push_back(h2);
    }

    TH2D* H2_wire_occupancy = new TH2D("wire_occupancy", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);

    /// --- Start analysis
    int nfile = 0;
    long unsigned int nevents =0;
    long unsigned int nelectrons =0;
    for (auto file : filenames) {
        nfile++;
        printf("> Open file %d/%d : %s\n", nfile, (int) filenames.size(), file.c_str());
        
        /// --- Initialise HIPO reader
        hipo::reader  reader(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        /// --- Define banks to be read
        hipo::bank  recBank(factory.getSchema("REC::Particle"));
        hipo::bank  adcBank(factory.getSchema("AHDC::adc"));


        /// --- Loop over events
        hipo::event event;
        long unsigned int nevents_per_file = 0;

        while( reader.next()){
            nevents++;
            nevents_per_file++;

            // display progress Bar
            if ((nevents_per_file % 1000 == 0) || ((int) nevents_per_file == reader.getEntries())) {
                progressBar(100.0*nevents_per_file/reader.getEntries());
            }
            // load bank content for this event
            reader.read(event);
            event.getStructure(recBank);
            
            for (int row = 0; row < recBank.getRows(); row++) {

                /// --- select trigger electrons
                if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0){ 

                    double px = recBank.getFloat("px", row) * Units::GeV;
                    double py = recBank.getFloat("py", row) * Units::GeV;
                    double pz = recBank.getFloat("pz", row) * Units::GeV;
                    double p = sqrt(px*px + py*py + pz*pz);

                    // physics kinematics
                    double scattered_beam_energy = sqrt(p*p + electron_mass*electron_mass);
                    double nu = beam_energy - scattered_beam_energy;
                    double theta = acos(pz/p) * Units::rad;
                    double Q2 = 4*beam_energy*scattered_beam_energy*pow(sin(theta/2),2);
                    double W2 = pow(deuteron_mass,2) + 2*deuteron_mass*nu - Q2;

                    double phi = atan2(py, px) * Units::rad;
                    if (phi < 0) phi += 2*M_PI;

                    H1_W2->Fill(W2);

                    if (W2 > W2_min && W2 < W2_max) {
                        nelectrons++;
                        /// --- loop over AHDC::adc
                        event.getStructure(adcBank);
                        for (int row = 0; row < adcBank.getRows(); row++) {
                            int layer = adcBank.getByte("layer", row);
                            int wire  = adcBank.getShort("component", row);

                            int num = AhdcUtils::layer2number(layer);

                            H2_corr_wire_phi[num-1]->Fill(phi / Units::deg, wire);
                            H2_wire_occupancy->Fill(wire, num);

                            //printf("%lf ", phi / Units::deg);
                            
                        } // end loop over adcBank rows
                        break; // we only select one electron per event
                    }

                } // end selection of trigger electrons

            } // end loop over rec particle

        } // end loop over events for this file
        printf("\033[34m nevents : %ld \033[0m \n", nevents);


    } // end loop over input files

    

    /// Store histograms
    TFile *f = new TFile(output.c_str(), "RECREATE");
    
    
    H1_W2->Write("W2");
    showCuts(H1_W2, W2_min, W2_max)->Write("W2_showing_elastic_limits");

    for (int i = 1; i <= 8; i++) {
        int layer = AhdcUtils::number2layer(i);
        std::string name = Form("wire_disposition_layer_%d", layer);
        H2_corr_wire_phi[i-1]->Write(name.c_str());
    }

    H2_wire_occupancy->Write("wire_occupancy");
    
    // ALL IN ONE
    {
        TCanvas* canvas = new TCanvas();
        canvas->Divide(3,3);

        // occupancy
        canvas->cd(1);
        H2_wire_occupancy->Draw("colz");
        
        // wire disposition
        for (int i = 1; i <= 8; i++) {
            canvas->cd(i+1);
            H2_corr_wire_phi[i-1]->Draw("colz");
        }

        canvas->Write("wire_disposition_all_in_one");

    }




    f->Close();
    printf("nb elastic electrons : %ld\n", nelectrons);
    printf("File created : %s\n", output.c_str());

    /// --- end of the program
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

TCanvas* showCuts(TH1D* h, double lim_inf, double lim_sup) {
    TCanvas* canvas = new TCanvas();
    canvas->cd();
    h->Draw();
    double ymax = h->GetMaximum();
    TLine* line_inf = new TLine(lim_inf, 0, lim_inf, ymax);
    TLine* line_sup = new TLine(lim_sup, 0, lim_sup, ymax);
    
    line_inf->SetLineColor(kRed);
    line_sup->SetLineColor(kRed);
    line_inf->SetLineWidth(2);
    line_sup->SetLineWidth(2);
    
    line_inf->Draw("same L");
    line_sup->Draw("same L");

    return canvas;
}