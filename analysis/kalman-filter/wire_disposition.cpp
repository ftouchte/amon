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
#include <map>

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
#include "TLatex.h"
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
TCanvas* xyPlotBadWires();
TCanvas* saveResidualLRBadWiresGr1(std::map<int, TH1D*> map_gr1);
void save_all_residual_LR_histos(std::vector<TH1D*> histos, TFile* file);

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
        h2->SetStats(false);
        H2_corr_wire_phi.push_back(h2);
    }

    TH2D* H2_wire_occupancy = new TH2D("wire_occupancy", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);
    H2_wire_occupancy->SetStats(false);

    /// --- Study bad wires
    std::map<int,TH1D*> bad_wires_map1;
    std::map<int,TH1D*> bad_wires_map2;
    std::map<int,TH1D*> bad_wires_map3;
    std::vector<int> bad_wires_vec1 = {24, 25, 75, 76, 77, 131, 132, 133};
    std::vector<int> bad_wires_vec2 = {34, 35, 88, 89, 90, 144, 145, 214};
    std::vector<int> bad_wires_vec3 = {47, 46, 102, 100, 226, 298, 157};
        // group 1
    for (int i = 0; i < (int) bad_wires_vec1.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires_vec1[i]);
        //int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        std::string name = Form("L%dW%d_bad_residual_LR", 10*sl+l, w);
        std::string title = Form("Num %d, L%dW%d, residual LR; residual LR (mm); count", bad_wires_vec1[i], 10*sl+l, w);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, -2.2, 2.2);
        
        bad_wires_map1[bad_wires_vec1[i]] = h;
    }
        // group 2
    for (int i = 0; i < (int) bad_wires_vec2.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires_vec2[i]);
        //int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        std::string name = Form("L%dW%d_bad_residual_LR", 10*sl+l, w);
        std::string title = Form("Num %d, L%dW%d, residual LR; residual LR (mm); count", bad_wires_vec2[i], 10*sl+l, w);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, -2.2, 2.2);
        
        bad_wires_map2[bad_wires_vec2[i]] = h;
    }
    // group 3
    for (int i = 0; i < (int) bad_wires_vec3.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires_vec3[i]);
        //int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        std::string name = Form("L%dW%d_bad_residual_LR", 10*sl+l, w);
        std::string title = Form("Num %d, L%dW%d, residual LR; residual LR (mm); count", bad_wires_vec3[i], 10*sl+l, w);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, -2.2, 2.2);
        
        bad_wires_map3[bad_wires_vec3[i]] = h;
    }

    
    /// ---
    std::vector<TH1D*> H1_all_residual_LR;
    for (int i = 0; i < 576; i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(i);
        int layer = ids[1];
        int wire = ids[2];
        std::string name = Form("L%dW%d_residual_LR", layer, wire);
        std::string title = Form("Num %d, L%dW%d, residual LR; residual LR (mm); count",i, layer, wire);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, -2.2, 2.2);
        H1_all_residual_LR.push_back(h);
    }

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
        hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
        hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));


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
                        for (int i = 0; i < adcBank.getRows(); i++) {
                            int layer = adcBank.getByte("layer", i);
                            int wire  = adcBank.getShort("component", i);
                            //int adc = adcBank.getInt("ADC", i);
                            //int wfType = adcBank.getInt("wfType", i);

                            //if (adc < 30) continue;
                            //if (wfType > 2) continue;

                            int num = AhdcUtils::layer2number(layer);

                            H2_corr_wire_phi[num-1]->Fill(phi / Units::deg, wire);
                            H2_wire_occupancy->Fill(wire, num);
                            
                        } // end loop over adcBank rows

                        /// --- find elastic tracks
                        event.getStructure(hitBank);
                        event.getStructure(trackBank);
                        int track_row = -1; 
                        for (int t = 0; t < trackBank.getRows(); t++) {
                            int nhits = trackBank.getInt("n_hits", t);
                            double track_px = trackBank.getFloat("px", row) * Units::MeV;
                            double track_py = trackBank.getFloat("py", row) * Units::MeV;
                            double track_phi = atan2(track_py, track_px) * Units::rad;
                            if (track_phi < 0) track_phi += 2*M_PI;
                            double delta_phi = fabs(fabs(phi-track_phi)-M_PI);
                            if (delta_phi < 20*Units::deg && nhits >= 7){ 
                                track_row = t;
                                break;
                            }
                        } // end loop over track

                        if (track_row < 0) continue;
                        int trackid = trackBank.getInt("trackid", track_row);
                        for (int i = 0; i < hitBank.getRows(); i++) {

                            if (trackid != hitBank.getInt("trackid", i)) continue;

                            int layer = 10*hitBank.getByte("superlayer", i) + hitBank.getByte("layer", i);
                            int wire = hitBank.getInt("wire", i);
                            double residual_LR = hitBank.getDouble("timeOverThreshold", i); // mm, contains residual LR
                            
                            int wireNum = AhdcUtils::slc2wire(1, layer, wire);

                            // all residual_LR
                            H1_all_residual_LR[wireNum]->Fill(residual_LR);

                            // bad wires gr 1
                            auto it1 = bad_wires_map1.find(wireNum);
                            if (it1 != bad_wires_map1.end()) {
                                it1->second->Fill(residual_LR);
                            }
                            // bad wires gr 2
                            auto it2 = bad_wires_map2.find(wireNum);
                            if (it2 != bad_wires_map2.end()) {
                                it2->second->Fill(residual_LR);
                            }
                            // bad wires gr 3
                            auto it3 = bad_wires_map3.find(wireNum);
                            if (it3 != bad_wires_map3.end()) {
                                it3->second->Fill(residual_LR);
                            }

                        } // end loop over hits

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

    xyPlotBadWires()->Write("xy_view_of_bad_wires");

    saveResidualLRBadWiresGr1(bad_wires_map1)->Write("residual_LR_bad_wires_gr_1");
    saveResidualLRBadWiresGr1(bad_wires_map2)->Write("residual_LR_bad_wires_gr_2");
    saveResidualLRBadWiresGr1(bad_wires_map3)->Write("residual_LR_bad_wires_gr_3");

    save_all_residual_LR_histos(H1_all_residual_LR, f);




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


TCanvas* saveResidualLRBadWiresGr1(std::map<int, TH1D*> map_gr1) {
    TCanvas* canvas = new TCanvas();
    canvas->Divide(3,3);
    int counter = 0;
    for (const auto& [wireNum, h] : map_gr1) {
        counter++;
        canvas->cd(counter);
        h->Draw();
    }
    return canvas;
}


#include "AhdcDetector.h"

TCanvas* xyPlotBadWires() {
    AhdcDetector* ahdc = new AhdcDetector();

    TGraph* gr_all_wires = new TGraph();
    gr_all_wires->SetTitle(";x (mm); y(mm)");
    gr_all_wires->SetMarkerStyle(4);
    gr_all_wires->SetMarkerColor(kBlue-9);

    for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
                    gr_all_wires->AddPoint(wire->x, wire->y);
                }
            }
        }
    }

    TGraph* gr_bad_wires = new TGraph();
    gr_bad_wires->SetTitle(";x (mm); y(mm)");
    gr_bad_wires->SetMarkerStyle(8);
    gr_bad_wires->SetMarkerColor(kRed);
    

    std::vector<int> bad_wires_id = {24, 25, 34, 35, 46, 47, 64, 65, 75, 76, 77, 88, 89, 90, 100, 102, 122, 124, 120, 131, 132, 133, 144, 145, 165, 167, 214, 226, 298};

    std::vector<TLatex*> bad_wires_labels;

    for (int i = 0; i < (int) bad_wires_id.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires_id[i]);
        int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        AhdcWire* wire = ahdc->GetSector(s-1)->GetSuperLayer(sl-1)->GetLayer(l-1)->GetWire(w-1);
        gr_bad_wires->AddPoint(wire->x, wire->y);

        TLatex* txt = new TLatex(wire->x, wire->y-0.3, TString::Format("%d, L%dW%d", bad_wires_id[i], ids[1], ids[2]).Data());
        txt->SetTextSize(0.01);
        txt->SetTextAlign(12);
        bad_wires_labels.push_back(txt);
    }

    TGraph* gr_fixed_wires = new TGraph();
    gr_fixed_wires->SetTitle(";x (mm); y(mm)");
    gr_fixed_wires->SetMarkerStyle(8);
    gr_fixed_wires->SetMarkerColor(kBlue);

    std::vector<int> fixed_wires_id = {64, 65, 122, 124, 120};

    for (int i = 0; i < (int) fixed_wires_id.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(fixed_wires_id[i]);
        int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        AhdcWire* wire = ahdc->GetSector(s-1)->GetSuperLayer(sl-1)->GetLayer(l-1)->GetWire(w-1);
        gr_fixed_wires->AddPoint(wire->x, wire->y);
    }


    // Rendering
    TCanvas* canvas = new TCanvas();
    canvas->cd();
    gr_all_wires->Draw("P");
    gr_bad_wires->Draw("P");
    gr_fixed_wires->Draw("P");

    // show text
    for (int i = 0; i < (int) bad_wires_labels.size(); i++) {
        bad_wires_labels[i]->Draw();
    }

    return canvas;

}

void save_all_residual_LR_histos(std::vector<TH1D*> histos, TFile* file) {
    // initialize canvas
    std::vector<TCanvas*> canvas;
    for (int i = 0; i < 48; i++) {
        TCanvas* c = new TCanvas(TString::Format("canvas_residual_LR_w%d-%d", 12*i, 12*(i+1)-1).Data(), "", 1500, 1600);
        c->Divide(3,4);
        canvas.push_back(c);
    }

    // draw histograms
    for (int i = 0; i < 576; i++) {
        int num = i / 12;
        int pad = i % 12 + 1;
        canvas[num]->cd(pad);
        histos[i]->Draw();
    }

    // save histos
    TDirectory * dir = file->mkdir("all_residuals_LR");
    dir->cd();
    for (int i = 0; i < 48; i++) {
        canvas[i]->Write(TString::Format("residual_LR_w%d-%d", 12*i, 12*(i+1)-1).Data());
        if (i == 0) {
            canvas[i]->Print("all_residual.pdf[");
        } 
        else if (i == 47) {
            canvas[i]->Print("all_residual.pdf]");
        } 
        else {
            canvas[i]->Print("all_residual.pdf");
        }
    }
    // go back to file
    file->cd();

}
