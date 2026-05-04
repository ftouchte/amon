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

#include "wire_disposition.h"
 
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
#include "TTree.h"

#include "futils.h"
#include "fOptions.h"
#include "Units.h"
#include "AhdcUtils.h"
#include "AhdcDetector.h"
#include "AhdcCCDB.h"


void progressBar(int state, int bar_length = 100);
TCanvas* showCuts(TH1D* h, double lim_inf, double lim_sup);
TCanvas* xyPlotBadWires(std::vector<int> & bad_wires, std::vector<int> & fixed_wires);
TCanvas* saveResidualLRBadWiresGroup(std::map<int, TH1D*> map_gr1, const char *name = nullptr);
void save_all_wires_histos(std::vector<TH1D*> histos, TFile* file, const char * name);
void plot_t0_distribution(TFile* file, std::vector<int> bad_wires);
void saveMappedHistos(std::map<int, TH1D*> map, TFile* file, const char * name);
std::map<int, TH1D*> createMapOfHistoWires(std::vector<int> bad_wires, const char * name, const char* title, int nbins, double xmin, double xmax);

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
    std::vector<int> all_bad_wires = {24, 25, 34, 35, 46, 47, 75, 76, 77, 88, 89, 90, 100, 102, 131, 132, 133, 140, 143, 144, 145, 146, 157, 165, 167, 196, 198, 214, 228, 369};
    std::vector<int> fixed_wires = {64, 65, 122, 124, 120};
    
    std::vector<int> bad_wires_vec1 = {24, 25, 75, 76, 77, 131, 132, 133};
    std::vector<int> bad_wires_vec2 = {34, 35, 88, 89, 90, 144, 145, 214};
    std::vector<int> bad_wires_vec3 = {47, 46, 102, 100, 226, 298, 157, 45, 158, 101};

    std::map<int,TH1D*> bad_wires_map1 = createMapOfHistoWires(bad_wires_vec1, "Gr1_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);
    std::map<int,TH1D*> bad_wires_map2 = createMapOfHistoWires(bad_wires_vec2, "Gr2_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);
    std::map<int,TH1D*> bad_wires_map3 = createMapOfHistoWires(bad_wires_vec3, "Gr3_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);

    std::map<int, TH1D*> leadingEdgeTimeBadWires_map = createMapOfHistoWires(all_bad_wires, "bad_leadingEdgeTime", "leadingEdgeTime; leadingEdgeTime (ns); count", 100, 200, 700);
    std::map<int,TH1D*>  residual_LRBadWires_map = createMapOfHistoWires(all_bad_wires, "bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);
    
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

    std::vector<TH1D*> H1_all_leadingEdgeTime;
    for (int i = 0; i < 576; i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(i);
        int layer = ids[1];
        int wire = ids[2];
        std::string name = Form("bad_L%dW%d_leadingEdgeTime", layer, wire);
        std::string title = Form("Num %d, L%dW%d, leadingEdgeTime; leadingEdgeTime (ns); count", i, layer, wire);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, 0, 700);
        H1_all_leadingEdgeTime.push_back(h);
    }

    /// --- Define TTree
    gInterpreter->GenerateDictionary("Hit;Track;std::vector<Hit>", "vector");
    TTree *tree = new TTree("tree", "tracks avec hits");

    Track track;

    tree->Branch("track", &track);

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
                            double track_pz = trackBank.getFloat("pz", row) * Units::MeV;
                            double track_p = sqrt(track_px*track_px + track_py*track_py + track_pz*track_pz);
                            double track_theta = acos(track_pz/track_p) * Units::rad;
                            double track_phi = atan2(track_py, track_px) * Units::rad;
                            if (track_phi < 0) track_phi += 2*M_PI;
                            double delta_phi = fabs(fabs(phi-track_phi)-M_PI);
                            if (delta_phi < 20*Units::deg && nhits >= 7){ 
                                track_row = t;
                                // fill TTree
                                track.px = px;
                                track.py = py;
                                track.pz = pz;
                                track.p = track_p;
                                track.theta = track_theta / Units::deg;
                                track.phi = track_phi / Units::deg;
                                track.nhits = trackBank.getInt("n_hits", t);

                                track.electron_p = p;
                                track.electron_px = px;
                                track.electron_py = py;
                                track.electron_pz = pz;
                                track.electron_theta = theta / Units::deg;
                                track.electron_phi = phi / Units::deg;
                                track.hits.clear();
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

                            {
                                auto it = residual_LRBadWires_map.find(wireNum);
                                if (it != residual_LRBadWires_map.end()) {
                                    it->second->Fill(residual_LR);
                                }
                            }

                            // leadingEdgeTime
                            int adcRow = hitBank.getShort("id", i) - 1;
                            double leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", adcRow);
                            H1_all_leadingEdgeTime[wireNum]->Fill(leadingEdgeTime);
                            {
                                auto it = leadingEdgeTimeBadWires_map.find(wireNum);
                                if (it != leadingEdgeTimeBadWires_map.end()) {
                                    it->second->Fill(leadingEdgeTime);
                                }
                            }

                            // Fill TTree
                            Hit h;
                            h.sector = 1;
                            h.superlayer = layer / 10;
                            h.layer = layer % 10;
                            h.component = wire;
                            h.wireNum = wireNum;
                            h.adc = hitBank.getInt("adc", i);
                            h.time = hitBank.getDouble("time", i);
                            h.residual_LR = residual_LR;
                            h.residual = hitBank.getDouble("residual", i);
                            track.hits.push_back(h);

                        } // end loop over hits

                        tree->Fill();

                        break; // we only select one electron per event
                    } // end cut on W2 to select elastics

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

    // Others
    xyPlotBadWires(all_bad_wires, fixed_wires)->Write("xy_view_of_bad_wires");

    saveResidualLRBadWiresGroup(bad_wires_map1, "residual_LR_bad_wires_gr_1.pdf")->Write("residual_LR_bad_wires_gr_1");
    saveResidualLRBadWiresGroup(bad_wires_map2, "residual_LR_bad_wires_gr_2.pdf")->Write("residual_LR_bad_wires_gr_2");
    saveResidualLRBadWiresGroup(bad_wires_map3, "residual_LR_bad_wires_gr_3.pdf")->Write("residual_LR_bad_wires_gr_3");

    save_all_wires_histos(H1_all_residual_LR, f, "all_wires_residual_LR");
    save_all_wires_histos(H1_all_leadingEdgeTime, f, "all_wires_leadingEdgeTime");

    saveMappedHistos(leadingEdgeTimeBadWires_map, f, "all_bad_wires_leadingEdgeTime");
    saveMappedHistos(residual_LRBadWires_map, f, "all_bad_wires_residual_LR");

    plot_t0_distribution(f, all_bad_wires);

    //tree->Write();


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

/**
 * @brief 
 * 
 * @param h histo 1D 
 * @param lim_inf where to put the first vertical line. Use -1e10 to exclude this limit
 * @param lim_sup where to put the second vertical line. Use +1e10 to exclude this limit
 * @return TCanvas* 
 */
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


TCanvas* saveResidualLRBadWiresGroup(std::map<int, TH1D*> map_gr1, const char *name) {
    TCanvas* canvas = new TCanvas(name, "", 1500, 1200);
    if (map_gr1.size() < 10) {
        canvas->Divide(3,3);
    } else {
        canvas->Divide(3,4);
    }
    
    int counter = 0;
    for (const auto& [wireNum, h] : map_gr1) {
        counter++;
        canvas->cd(counter);
        h->Draw();
    }
    if (name != nullptr)
        canvas->Print(TString::Format("./output/%s.pdf", name).Data());
    return canvas;
}

/**
 * @brief 
 * 
 * @param bad_wires List of unique wire id (numerotation starting at 0)
 * @param fixed_wires List of unique wire id (numerotation starting at 0)
 * @return TCanvas* conating the plot ready to be stored or printed
 */
TCanvas* xyPlotBadWires(std::vector<int> & bad_wires, std::vector<int> & fixed_wires) {
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

    std::vector<TLatex*> bad_wires_labels;

    for (int i = 0; i < (int) bad_wires.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires[i]);
        int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        AhdcWire* wire = ahdc->GetSector(s-1)->GetSuperLayer(sl-1)->GetLayer(l-1)->GetWire(w-1);
        gr_bad_wires->AddPoint(wire->x, wire->y);

        TLatex* txt = new TLatex(wire->x, wire->y-0.3, TString::Format("%d, L%dW%d", bad_wires[i], ids[1], ids[2]).Data());
        txt->SetTextSize(0.01);
        txt->SetTextAlign(12);
        bad_wires_labels.push_back(txt);
    }

    TGraph* gr_fixed_wires = new TGraph();
    gr_fixed_wires->SetTitle(";x (mm); y(mm)");
    gr_fixed_wires->SetMarkerStyle(8);
    gr_fixed_wires->SetMarkerColor(kBlue);


    for (int i = 0; i < (int) fixed_wires.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(fixed_wires[i]);
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

/**
 * @brief 
 * 
 * @param histos vector of 576 histos
 * @param file pointer of TFile*
 * @param name name the quantity to plot be plotted
 */
void save_all_wires_histos(std::vector<TH1D*> histos, TFile* file, const char * name) {
    // initialize canvas
    std::vector<TCanvas*> canvas;
    for (int i = 0; i < 48; i++) {
        TCanvas* c = new TCanvas(TString::Format("canvas_%s_w%d-%d", name, 12*i, 12*(i+1)-1).Data(), "", 1500, 1600);
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
    TDirectory * dir = file->mkdir(TString::Format("%s[", name).Data());
    dir->cd();

    canvas[0]->Print(TString::Format("./output/%s.pdf[", name).Data()); // open file
    for (int i = 0; i < 48; i++) {
        canvas[i]->Write(TString::Format("%s_w%d-%d", name, 12*i, 12*(i+1)-1).Data());
        canvas[i]->Print(TString::Format("./output/%s.pdf", name).Data());
    }
    canvas[47]->Print(TString::Format("./output/%s.pdf]", name).Data()); // close file

    // go back to file
    file->cd();
}

/**
 * @brief 
 * 
 * @param map of histograms 
 * @param file pointer to the TFile*
 * @param name name of the histogram 1D
 */
void saveMappedHistos(std::map<int, TH1D*> map, TFile* file, const char * name) {

    // for printing
    int size = map.size();
    int nb_canvas = (size/12) + (size % 12 == 0 ? 0 : 1);
    std::vector<TCanvas*> canvas;
    for (int i = 0; i < nb_canvas; i++) {
        TCanvas* c = new TCanvas(TString::Format("canvas_%s_w%d-%d", name, 12*i, 12*(i+1)-1).Data(), "", 1500, 1600);
        c->Divide(3,4);
        canvas.push_back(c);
    }

    // save in root file
    TDirectory * dir = file->mkdir(name);
    dir->cd();

    int counter = 0;
    for (const auto& [wireNum, h] : map) {
        // writing
        h->Write();

        // printing
        int num = counter / 12;
        int pad = counter % 12 + 1;
        canvas[num]->cd(pad);
        h->Draw();
        counter++;
    }

    // go back to file
    file->cd();

    canvas[0]->Print(TString::Format("./output/%s.pdf[", name).Data()); // open file
    for (int i = 0; i < nb_canvas; i++) {
        canvas[i]->Print(TString::Format("./output/%s.pdf", name).Data());
    }
    canvas[nb_canvas-1]->Print(TString::Format("./output/%s.pdf]", name).Data()); // close file
}

void plot_t0_distribution(TFile* file, std::vector<int> bad_wires) {

    AhdcCCDB* ccdb = new AhdcCCDB("mysql://clas12reader@clasdb.jlab.org/clas12", 22712, "default", "no");

    TH1D* h = new TH1D("t0", "wire t0 distribution; t0 (ns); count", 100, 150, 380);
    
    TGraph* gr = new TGraph(576);
    TGraph* gr2 = new TGraph(bad_wires.size());
    // gr->SetTitle("t0 versus wire; wire; t0 (ns)");
    gr->SetMarkerStyle(20);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kRed);

    int counter = 0;
    for (int i = 0; i < 576; i++) {
        double t0 = ccdb->get_t0(i).t0;
        //double dt0 = ccdb->get_t0(i).dt0;
        h->Fill(t0);
        gr->SetPoint(i, i, t0);
        if (t0 > 220) {
            counter++;
            std::vector<int> ids = AhdcUtils::wire2slc(i);
            printf("%3d  |  %3d    %2d   %2d   %lf\n", i+1, i, ids[1], ids[2], t0);
        }
    }

    for (int i = 0; i < (int) bad_wires.size(); i++) {
        double t0 = ccdb->get_t0(bad_wires[i]).t0;
        //double dt0 = ccdb->get_t0(i).dt0;
        h->Fill(t0);
        gr2->SetPoint(i, bad_wires[i], t0);
    }

    printf("(nwires   :   %d)\n", counter);

    TDirectory * dir = file->mkdir("t0_distribution");
    dir->cd();
    
    h->Write("t0");
    //gr->Write("t0_versus_wire");

    TCanvas* c = new TCanvas();
    c->cd();
    gr->Draw("AP");
    gr2->Draw("same p");
    c->Write("t0_versus_wire");
    

    // go back to file
    file->cd();
}



/**
 * @brief Create a Map Of Bad Wires object
 * 
 * @param bad_wires List of unique wire id (numerotation starting at 0)
 * @return std::map<int, TH1D*> 
 */
std::map<int, TH1D*> createMapOfHistoWires(std::vector<int> bad_wires, const char * _name, const char* _title, int nbins, double xmin, double xmax) {

    std::map<int, TH1D*> map;

    for (int i = 0; i < (int) bad_wires.size(); i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(bad_wires[i]);
        //int s = ids[0];
        int sl = ids[1] / 10;
        int l = ids[1] % 10;
        int w = ids[2];
        std::string name = Form("L%dW%d_%s", 10*sl+l, w, _name);
        std::string title = Form("Num %d, L%dW%d, %s", bad_wires[i], 10*sl+l, w, _title);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
        
        map[bad_wires[i]] = h;
    }

    return map;

}

//# From the wire swap study we found that bad wires have a very large t0. And for those that we managed to swap, they had a good t0. Assuming that the fit was not good, we set their t0 to the mean of the other wires
//# Sector Layer Component t0 dt0 extra1 extra2 chi2ndf
// t_{0} versus wire , run 22712, timestamp : May 04, 2026