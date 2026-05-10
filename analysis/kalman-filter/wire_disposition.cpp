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
#include <thread>

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
std::vector<TH2D*> generate_corr_wire_phi(const char* tag);
TCanvas* save_wire_versus_phi_plots (std::vector<TH2D*> H2_corr_wire_phi, TH2D* H2_wire_occupancy);

/// --- Constants
//const double beam_energy =  10676.6 * Units::MeV;
const double beam_energy = 2.23951 * Units::GeV; // incident energy if the electron, GeV
const double electron_mass = 0.511e-3 * Units::GeV; // energy mass of electron, GeV
const double proton_mass = 938.272e-3 * Units::GeV; // energy mass of proton, GeV
const double helium_mass = 3.73 * Units::GeV; // energy mass of Helium-4, GeV
const double deuteron_mass = 1.875 * Units::GeV; // energy mass of Deuterium, GeV

int main(int argc, char const *argv[]) {

    /// --- start timer
    auto start = std::chrono::high_resolution_clock::now();

    /// --- Load options
    fOptions OPT({"-i", "-o", "-v", "-simu"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    SingleThreadRun t;

    /// --- Parallel computing
    int num_threads = 10;
    std::vector<std::thread> threads(num_threads-1);
    std::vector<SingleThreadRun> others_t(num_threads-1);

    int num_files = filenames.size();
    int num_files_per_threads = num_files / num_threads + (num_files % num_threads == 0 ? 0 : 1);
    auto iterator = filenames.begin();
    
    for (int i = 0; i < num_threads-1; i++) {
        auto iterator_end = std::next(iterator, num_files_per_threads);
        //std::vector<std::string> _filenames(iterator, iterator_end);

        threads[i] = std::thread(&SingleThreadRun::run, &others_t[i], std::vector<std::string>(iterator, iterator_end));

        iterator = iterator_end;
    }

    // the current thread processes the rest of the files
    t.run(std::vector<std::string>(iterator, filenames.end()));
    
    // wait for the others threads to join
    for (auto& thr : threads) {
        thr.join();
    }

    // now merge execution
    for (auto & ti : others_t) {
        t.merge(ti);
    }

    

    /// --- Output
    TFile *f = new TFile(output.c_str(), "RECREATE");
    
    
    t.H1_W2->Write("W2");
    showCuts(t.H1_W2, t.W2_min, t.W2_max)->Write("W2_showing_elastic_limits");

    t.H1_track_vz->Write("track_vz");
    showCuts(t.H1_track_vz, t.vzInf_min, t.vzInf_max)->Write("track_vzInf_limits");
    showCuts(t.H1_track_vz, t.vzSup_min, t.vzSup_max)->Write("track_vzSup_limits");

    save_wire_versus_phi_plots(t.H2_corr_wire_phi, t.H2_wire_occupancy)->Write("wire_disposition_all_in_one");
    save_wire_versus_phi_plots(t.H2_corr_wire_phi_vzInf, t.H2_wire_occupancy_vzInf)->Write("wire_disposition_all_in_one_vzInf");
    save_wire_versus_phi_plots(t.H2_corr_wire_phi_vzSup, t.H2_wire_occupancy_vzSup)->Write("wire_disposition_all_in_one_vzSup");

    // Others
    xyPlotBadWires(t.all_bad_wires, t.fixed_wires)->Write("xy_view_of_bad_wires");

    saveResidualLRBadWiresGroup(t.bad_wires_map1, "residual_LR_bad_wires_gr_1.pdf")->Write("residual_LR_bad_wires_gr_1");
    saveResidualLRBadWiresGroup(t.bad_wires_map2, "residual_LR_bad_wires_gr_2.pdf")->Write("residual_LR_bad_wires_gr_2");
    saveResidualLRBadWiresGroup(t.bad_wires_map3, "residual_LR_bad_wires_gr_3.pdf")->Write("residual_LR_bad_wires_gr_3");

    save_all_wires_histos(t.H1_all_residual_LR, f, "all_wires_residual_LR");
    save_all_wires_histos(t.H1_all_leadingEdgeTime, f, "all_wires_leadingEdgeTime");

    saveMappedHistos(t.leadingEdgeTimeBadWires_map, f, "all_bad_wires_leadingEdgeTime");
    saveMappedHistos(t.residual_LRBadWires_map, f, "all_bad_wires_residual_LR");

    plot_t0_distribution(f, t.all_bad_wires);

    f->Close();
    printf("nb elastic electrons : %ld\n", t.nelectrons);
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


std::vector<TH2D*> generate_corr_wire_phi(const char* tag) {
    std::vector<TH2D*> H2_corr_wire_phi;
    for (int i = 1; i <= 8; i++) {
        int layer = AhdcUtils::number2layer(i);
        int nbWires = AhdcUtils::layerNbWires(layer);
        std::string name = Form("wire_disposition_layer_%d%s", layer, tag);
        std::string title = Form("Wire occupation versus Electron phi (layer %d); #phi_{e} (deg); wire", layer);
        TH2D* h2 = new TH2D(name.c_str(), title.c_str(), nbWires, 0, 360, nbWires, 1, nbWires+1);
        h2->SetStats(false);
        H2_corr_wire_phi.push_back(h2);
    }
    return H2_corr_wire_phi;
}

TCanvas* save_wire_versus_phi_plots (std::vector<TH2D*> H2_corr_wire_phi, TH2D* H2_wire_occupancy) {
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

    return canvas;
    //canvas->Write("wire_disposition_all_in_one");
}


//# From the wire swap study we found that bad wires have a very large t0. And for those that we managed to swap, they had a good t0. Assuming that the fit was not good, we set their t0 to the mean of the other wires
//# Sector Layer Component t0 dt0 extra1 extra2 chi2ndf
// t_{0} versus wire , run 22712, timestamp : May 04, 2026

SingleThreadRun::SingleThreadRun () {

    H1_W2 = new TH1D("W2", "W^{2}; W^{2} (GeV^{2}); count", 100, 3.2, 6);
    H1_track_vz = new TH1D("track_vz", "track vz; vz (cm); count", 100, -24, 20);

    H2_corr_wire_phi = generate_corr_wire_phi("");
    H2_corr_wire_phi_vzInf = generate_corr_wire_phi("_vzInf");
    H2_corr_wire_phi_vzSup = generate_corr_wire_phi("_vzSup");

    H2_wire_occupancy = new TH2D("wire_occupancy", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);
    H2_wire_occupancy->SetStats(false);
    H2_wire_occupancy_vzInf = new TH2D("wire_occupancy_vzInf", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);
    H2_wire_occupancy_vzInf->SetStats(false);
    H2_wire_occupancy_vzSup = new TH2D("wire_occupancy_vzSup", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);
    H2_wire_occupancy_vzSup->SetStats(false);

    bad_wires_map1 = createMapOfHistoWires(bad_wires_vec1, "Gr1_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);
    bad_wires_map2 = createMapOfHistoWires(bad_wires_vec2, "Gr2_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);
    bad_wires_map3 = createMapOfHistoWires(bad_wires_vec3, "Gr3_bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);

    leadingEdgeTimeBadWires_map = createMapOfHistoWires(all_bad_wires, "bad_leadingEdgeTime", "leadingEdgeTime; leadingEdgeTime (ns); count", 100, 200, 700);
    residual_LRBadWires_map = createMapOfHistoWires(all_bad_wires, "bad_residual_LR", "residul LR; residual LR (mm); count", 100, -2.2, 2.2);

    for (int i = 0; i < 576; i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(i);
        int layer = ids[1];
        int wire = ids[2];
        std::string name = Form("L%dW%d_residual_LR", layer, wire);
        std::string title = Form("Num %d, L%dW%d, residual LR; residual LR (mm); count",i, layer, wire);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, -2.2, 2.2);
        H1_all_residual_LR.push_back(h);
    }

    for (int i = 0; i < 576; i++) {
        std::vector<int> ids = AhdcUtils::wire2slc(i);
        int layer = ids[1];
        int wire = ids[2];
        std::string name = Form("bad_L%dW%d_leadingEdgeTime", layer, wire);
        std::string title = Form("Num %d, L%dW%d, leadingEdgeTime; leadingEdgeTime (ns); count", i, layer, wire);
        TH1D* h = new TH1D(name.c_str(), title.c_str(), 100, 0, 700);
        H1_all_leadingEdgeTime.push_back(h);
    }

    /// --- Stats
    nfiles = 0;
    nevents = 0;
    nelectrons = 0;

}

SingleThreadRun::~SingleThreadRun () {
    delete H1_W2;
    delete H1_track_vz;
    delete H2_wire_occupancy;
    delete H2_wire_occupancy_vzInf;
    delete H2_wire_occupancy_vzSup;

    delete_histo_vector(H2_corr_wire_phi);
    delete_histo_vector(H2_corr_wire_phi_vzInf);
    delete_histo_vector(H2_corr_wire_phi_vzSup);

    delete_map(bad_wires_map1);
    delete_map(bad_wires_map2);
    delete_map(bad_wires_map3);
    delete_map(leadingEdgeTimeBadWires_map);
    delete_map(residual_LRBadWires_map);

    delete_histo_vector(H1_all_residual_LR);
    delete_histo_vector(H1_all_leadingEdgeTime);
    
}

void SingleThreadRun::add_map_histo(std::map<int,TH1D*> & map_dest, const std::map<int,TH1D*> & map_src) {
    for (const auto& [wireNum, h] : map_src) {
        
        auto it = map_dest.find(wireNum);

        if (it != map_dest.end()) {
            it->second->Add(h);
        }
    }
}

void SingleThreadRun::merge(const SingleThreadRun & t) {
    this->H1_W2->Add(t.H1_W2);
    this->H1_track_vz->Add(t.H1_track_vz);
    this->H2_wire_occupancy->Add(t.H2_wire_occupancy);
    this->H2_wire_occupancy_vzInf->Add(t.H2_wire_occupancy_vzInf);
    this->H2_wire_occupancy_vzSup->Add(t.H2_wire_occupancy_vzSup);

    add_map_histo(this->bad_wires_map1, t.bad_wires_map1);
    add_map_histo(this->bad_wires_map2, t.bad_wires_map2);
    add_map_histo(this->bad_wires_map3, t.bad_wires_map3);

    add_map_histo(this->leadingEdgeTimeBadWires_map, t.leadingEdgeTimeBadWires_map);
    add_map_histo(this->residual_LRBadWires_map, t.residual_LRBadWires_map);

    add_histo_vector(this->H2_corr_wire_phi, t.H2_corr_wire_phi);
    add_histo_vector(this->H2_corr_wire_phi_vzInf, t.H2_corr_wire_phi_vzInf);
    add_histo_vector(this->H2_corr_wire_phi_vzSup, t.H2_corr_wire_phi_vzSup);

    add_histo_vector(this->H1_all_residual_LR, t.H1_all_residual_LR);
    add_histo_vector(this->H1_all_leadingEdgeTime, t.H1_all_leadingEdgeTime);

    /// --- Stats
    this->nfiles += t.nfiles;
    this->nelectrons += t.nelectrons;
    this->nevents += t.nevents;

}

void SingleThreadRun::run(std::vector<std::string> filenames) {
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

                        } // end loop over hits

                        double vz = trackBank.getFloat("z", track_row) * Units::mm;
                        H1_track_vz->Fill(vz);
                        for (int i = 0; i < adcBank.getRows(); i++) {
                            int layer = adcBank.getByte("layer", i);
                            int wire  = adcBank.getShort("component", i);
                            //int adc = adcBank.getInt("ADC", i);
                            //int wfType = adcBank.getInt("wfType", i);

                            //if (adc < 30) continue;
                            //if (wfType > 2) continue;

                            int num = AhdcUtils::layer2number(layer);

                            // vzInf
                            if (vz > vzInf_min && vz < vzInf_max) {
                                H2_corr_wire_phi_vzInf[num-1]->Fill(phi / Units::deg, wire);
                                H2_wire_occupancy_vzInf->Fill(wire, num);
                            }
                            // vzSup
                            if (vz > vzSup_min && vz < vzSup_max) {
                                H2_corr_wire_phi_vzSup[num-1]->Fill(phi / Units::deg, wire);
                                H2_wire_occupancy_vzSup->Fill(wire, num);
                            }
                        } // end loop over adc

                        break; // we only select one electron per event
                    } // end cut on W2 to select elastics

                } // end selection of trigger electrons

            } // end loop over rec particle

        } // end loop over events for this file
        printf("\033[34m nevents : %ld \033[0m \n", nevents);


    } // end loop over input files
}

