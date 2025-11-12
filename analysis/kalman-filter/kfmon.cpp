/***********************************************
 * Kalman Filter monitoring for AHDC
 *
 *
 * @author Felix Touchte Codjo
 * @date November 7, 2025
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <tuple>
#include <map>

#include <vector>
#include <string>
#include <chrono>

#include "reader.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "Math/PdfFuncMathCore.h"
#include "TApplication.h"
#include "TText.h"
#include "THStack.h"
#include "TMultiGraph.h"

#include "AhdcCCDB.h"
#include "futils.h"
#include "AhdcDetector.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
    const char * output = "./output/kfmon-v20.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v20.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
    hipo::bank  kfmonBank(factory.getSchema("AHDC::kftrack:mon"));
    hipo::bank  track0Bank(factory.getSchema("AHDC::track"));
    hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
    hipo::bank  mcBank(factory.getSchema("MC::Particle"));
    hipo::bank  trueBank(factory.getSchema("MC::True"));
    hipo::event event;
    long unsigned int nevents =0;
    long unsigned int nMCtracks =0;
    long unsigned int nKFtracks =0;
   
    // AHDC view
    AhdcDetector * ahdc = new AhdcDetector();
    TGraph * gr_ahdcView = new TGraph(577);
    gr_ahdcView->SetPoint(0, 0, 0);
    double zpos = 0;
    int nwires = 0;
    for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
                    nwires++;
                    // view at z = 0
                    zpos = 0;
                    wire->set_z(zpos);
                    gr_ahdcView->SetPoint(nwires, wire->x, wire->y);
                }
            }
        }
    }
    gr_ahdcView->SetTitle(TString::Format("AHDC XY view at z = %.2lf; x(mm); y(mm)", zpos).Data());
    gr_ahdcView->SetMarkerStyle(4);
    gr_ahdcView->SetMarkerColorAlpha(kBlue-10, 0.3);
    TCanvas * c0 = new TCanvas("c0", "c0 ahdcView", 1000, 1000);
    c0->cd();
    gr_ahdcView->Draw("AP");
    TFile *f = new TFile(output, "RECREATE");
    c0->Write("ahdcView_empty");
    // color palette
    //std::vector<int> ColorPalette= {};

    // Loop over events
    while( reader.next()){
        if (nevents > 10) break;
        nevents++;
        // Progress Bar
        //if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
        //    progressBar(100.0*nevents/reader.getEntries());
        //}
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(trackBank);
        event.getStructure(kfmonBank);
        event.getStructure(track0Bank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);
        event.getStructure(trueBank);

        nMCtracks += mcBank.getRows();
        nKFtracks += trackBank.getRows();

        // loop over track
        for (int trackRow = 0; trackRow < trackBank.getRows(); trackRow++) {
            ///////////////////////////////
            // Display : states/positions
            // ////////////////////////////
            TCanvas * cn = new TCanvas(TString::Format("c%d_%d", (int) nevents, trackRow+1).Data(),TString::Format("c%d_%d ahdcView", (int) nevents, trackRow+1).Data(), 1600, 900);
            cn->Divide(2,1);
            // XY view
            cn->cd(1);
            gr_ahdcView->Draw("AP"); // AHDC view
            int trackid = trackBank.getInt("trackid", trackRow);
                // Loop over hits
            TGraph * gr_hits = new TGraph();
            gr_hits->SetMarkerColor(2);
            gr_hits->SetMarkerStyle(8);
            for (int i = 0; i < hitBank.getRows(); i++) {
                if (hitBank.getInt("trackid", i) != trackid) continue;
                int sl = hitBank.get("superlayer", i)-1;
                int l = hitBank.get("layer", i)-1;
                int w = hitBank.get("wire", i)-1;
				AhdcWire* wire = ahdc->GetSector(0)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
                wire->set_z(zpos);
                gr_hits->AddPoint(wire->x, wire->y);
            }
            gr_hits->Draw("P");
            // ZR view
                // define the frame
            TGraph * gr_zrView = new TGraph();
            gr_zrView->SetTitle("ZR view; z (mm); r (mm)");
            gr_zrView->AddPoint(-150,-5);
            gr_zrView->AddPoint(-150,70);
            gr_zrView->AddPoint(150,70);
            gr_zrView->AddPoint(150,-5);
            cn->cd(2);
            gr_zrView->Draw("AP");
            //gr_zrView->SetLineColor(kBlue-10);
                // draw nivel line
            for (auto radius : std::vector<double>({0, 3.0, 3.06, 32, 38, 38+4, 48, 48+4, 58, 58+4, 68})) {
                TGraph * gr = new TGraph();
                gr->AddPoint(-150,radius);
                gr->AddPoint(150,radius);
                gr->SetLineColor(kBlue-10);
                gr->SetLineStyle(2);
                //gr->SetLineWidth(2);
                gr->Draw("same l");
            }
            // key : Niter, orientation, status --> TGraph* to draw markers
            std::map<std::tuple<int, int, int>, std::pair<TGraph*,TGraph*>> graphMap;
            // key : Niter to draw lines
            std::map<int, std::pair<TGraph*,TGraph*>> graphMapIter;
            // Variances of x, y, z, px, py, pz around the beam axis (indicator ==0) 
            //std::vector<TGraph*> graphVariance = {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()};
            std::vector<std::vector<TGraph*>> graphVariance = {
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // prediction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // correction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}  // both of them
            };
            std::vector<std::vector<TGraph*>> graphValue = {
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // prediction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // correction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()},  // both of them
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}  // mc value
            };
            for (int i = 0; i < kfmonBank.getRows(); i++) {
                if (kfmonBank.getInt("trackid", i) != trackid) continue;
                int Niter = kfmonBank.getShort("Niter", i);
                int orientation = kfmonBank.getShort("orientation", i);
                int indicator = kfmonBank.getShort("indicator", i);
                int status = kfmonBank.getShort("status", i);
                double x = kfmonBank.getFloat("x", i);
                double y = kfmonBank.getFloat("y", i);
                double z = kfmonBank.getFloat("z", i);
                double px = kfmonBank.getFloat("px", i);
                double py = kfmonBank.getFloat("py", i);
                double pz = kfmonBank.getFloat("pz", i);
                double var_x = kfmonBank.getFloat("var_x", i);
                double var_y = kfmonBank.getFloat("var_y", i);
                double var_z = kfmonBank.getFloat("var_z", i);
                double var_px = kfmonBank.getFloat("var_px", i);
                double var_py = kfmonBank.getFloat("var_py", i);
                double var_pz = kfmonBank.getFloat("var_pz", i);
                
                // it : iterator, is an object that points to a std::pair containing 
                // it->first is the key, here the std::tuple
                // it->second is the value, here the TGraph*
                // inserted is a boolean, is false if the element already exist and true if it is a new element
                // if false, "it" corresponds to the existing element
                
                if (Niter >= 10) continue;
                { // Niter, orientation, status
                    auto [it, inserted] = graphMap.insert({std::make_tuple(Niter, orientation, status), std::make_pair(new TGraph(), new TGraph())});
                    (it->second).first->AddPoint(x,y); // xy
                    (it->second).second->AddPoint(z,sqrt(x*x+y*y)); //zr
                }
                { // Niter only
                    auto [it, inserted] = graphMapIter.insert({Niter, std::make_pair(new TGraph(), new TGraph())});
                    (it->second).first->AddPoint(x,y); // xy
                    (it->second).second->AddPoint(z,sqrt(x*x+y*y)); //zr
                }
                { // Variance graph
                    if (indicator == 0) { // beam level
                       graphVariance[2][0]->AddPoint(Niter, var_x); 
                       graphVariance[2][1]->AddPoint(Niter, var_y); 
                       graphVariance[2][2]->AddPoint(Niter, var_z); 
                       graphVariance[2][3]->AddPoint(Niter, var_px); 
                       graphVariance[2][4]->AddPoint(Niter, var_py); 
                       graphVariance[2][5]->AddPoint(Niter, var_pz);
                       if (status == 0) { // prediction
                           graphVariance[0][0]->AddPoint(Niter, var_x); 
                           graphVariance[0][1]->AddPoint(Niter, var_y); 
                           graphVariance[0][2]->AddPoint(Niter, var_z); 
                           graphVariance[0][3]->AddPoint(Niter, var_px); 
                           graphVariance[0][4]->AddPoint(Niter, var_py); 
                           graphVariance[0][5]->AddPoint(Niter, var_pz);
                       } 
                       else { // correction, status == 1
                           graphVariance[1][0]->AddPoint(Niter, var_x); 
                           graphVariance[1][1]->AddPoint(Niter, var_y); 
                           graphVariance[1][2]->AddPoint(Niter, var_z); 
                           graphVariance[1][3]->AddPoint(Niter, var_px); 
                           graphVariance[1][4]->AddPoint(Niter, var_py); 
                           graphVariance[1][5]->AddPoint(Niter, var_pz);
                       }
                    }
                }
                { // Value graph
                    if (indicator == 0) { // beam level
                       graphValue[2][0]->AddPoint(Niter, x); 
                       graphValue[2][1]->AddPoint(Niter, y); 
                       graphValue[2][2]->AddPoint(Niter, z); 
                       graphValue[2][3]->AddPoint(Niter, px); 
                       graphValue[2][4]->AddPoint(Niter, py); 
                       graphValue[2][5]->AddPoint(Niter, pz);
                       if (status == 0) { // prediction
                           graphValue[0][0]->AddPoint(Niter, x); 
                           graphValue[0][1]->AddPoint(Niter, y); 
                           graphValue[0][2]->AddPoint(Niter, z); 
                           graphValue[0][3]->AddPoint(Niter, px); 
                           graphValue[0][4]->AddPoint(Niter, py); 
                           graphValue[0][5]->AddPoint(Niter, pz);
                       } 
                       else { // correction, status == 1
                           graphValue[1][0]->AddPoint(Niter, x); 
                           graphValue[1][1]->AddPoint(Niter, y); 
                           graphValue[1][2]->AddPoint(Niter, z); 
                           graphValue[1][3]->AddPoint(Niter, px); 
                           graphValue[1][4]->AddPoint(Niter, py); 
                           graphValue[1][5]->AddPoint(Niter, pz);
                       }
                       graphValue[3][0]->AddPoint(Niter, mcBank.get("vx", 0)*10); 
                       graphValue[3][1]->AddPoint(Niter, mcBank.get("vy", 0)*10); 
                       graphValue[3][2]->AddPoint(Niter, mcBank.get("vz", 0)*10); 
                       graphValue[3][3]->AddPoint(Niter, mcBank.get("px", 0)*1000); 
                       graphValue[3][4]->AddPoint(Niter, mcBank.get("py", 0)*1000); 
                       graphValue[3][5]->AddPoint(Niter, mcBank.get("pz", 0)*1000);
                    }
                }
            }
            TLegend* legend = new TLegend(0.1,0.1,0.35,0.9);
            // Line style 
            for (auto& [k, g] : graphMapIter) {
                // Set Color
                // they should be the same as the marker
                int Niter = k;
                if (Niter <= 8) {
                    g.first->SetLineColor(Niter+1); 
                    g.second->SetLineColor(Niter+1);
                }
                else { // Niter > 9
                       // we ship the kColor 10 as it is transparente
                    g.first->SetLineColor(Niter+2);
                    g.second->SetLineColor(Niter+2);
                }
                // Draw line
                cn->cd(1);
                g.first->Draw("same L");
                cn->cd(2);
                g.second->Draw("same L");
            }
            // Marker style
            for (auto& [k, g] : graphMap) {
                auto [Niter, orientation, status] = k; // key
                // Set Color
                if (Niter <= 8) {
                    g.first->SetLineColor(Niter+1);
                    g.second->SetLineColor(Niter+1);
                    g.first->SetMarkerColor(Niter+1);
                    g.second->SetMarkerColor(Niter+1);
                }
                else { // Niter > 9
                       // we ship the kColor 10 as it is transparente
                    g.first->SetLineColor(Niter+2);
                    g.second->SetLineColor(Niter+2);
                    g.first->SetMarkerColor(Niter+2);
                    g.second->SetMarkerColor(Niter+2);
                }
                // Set Marker style
                g.first->SetMarkerSize(2);
                g.second->SetMarkerSize(2);
                if ((orientation % 2 == 0) && (status == 0)) {
                    // forward or postfit propagation and prediction
                    g.first->SetMarkerStyle(25);
                    g.second->SetMarkerStyle(25);
                    legend->AddEntry(g.first, TString::Format("Niter : %d, %s, %s", Niter, "forward", "prediction").Data());
                }
                else if ((orientation % 2 == 0) && (status == 1)) {
                    // forward or postfit propagation and correction
                    g.first->SetMarkerStyle(36);
                    g.second->SetMarkerStyle(36);
                    legend->AddEntry(g.first, TString::Format("Niter : %d, %s, %s", Niter, "forward", "correction").Data());
                }
                else if ((orientation % 2 == 1) && (status == 0)) {
                    // backward propagation and prediction
                    g.first->SetMarkerStyle(24);
                    g.second->SetMarkerStyle(24);
                    legend->AddEntry(g.first, TString::Format("Niter : %d, %s, %s", Niter, "backward", "prediction").Data());
                }
                else { //if ((orientation % 2 == 0) && (status == 1)) 
                    // backward propagation and correction
                    g.first->SetMarkerStyle(38);
                    g.second->SetMarkerStyle(38);
                    legend->AddEntry(g.first, TString::Format("Niter : %d, %s, %s", Niter, "backward", "correction").Data());
                }
                // XY view
                cn->cd(1);
                g.first->Draw("same P");
                // ZR view
                cn->cd(2);
                g.second->Draw("same P");
            }
            // Draw true info
            std::pair<TGraph*, TGraph*> gr_true(new TGraph(), new TGraph());
            { // MC::Particle
              // convert distance en mm
                double x = mcBank.get("vx", 0)*10;
                double y = mcBank.get("vy", 0)*10;
                double z = mcBank.get("vz", 0)*10;
                double r = sqrt(x*x+y*y);
                gr_true.first->AddPoint(x,y);
                gr_true.second->AddPoint(z,r);
            }
            for (int i = 0; i < trueBank.getRows(); i++) {
                int pid = trueBank.get("pid", i);
                if (pid != 1000010020) continue;
                // only look at deuteron
                double x = trueBank.get("avgX", i);
                double y = trueBank.get("avgY", i);
                double z = trueBank.get("avgZ", i);
                double r = sqrt(x*x+y*y);
                gr_true.first->AddPoint(x,y);
                gr_true.second->AddPoint(z,r);
            }
            gr_true.first->SetMarkerStyle(106);
            gr_true.first->SetMarkerSize(2);
            gr_true.first->SetMarkerColor(kRed);
            gr_true.first->SetLineColor(kRed);
            gr_true.first->SetLineStyle(1);
            gr_true.first->SetLineWidth(2);
            gr_true.second->SetMarkerStyle(106);
            gr_true.second->SetMarkerSize(2);
            gr_true.second->SetMarkerColor(kRed);
            gr_true.second->SetLineColor(kRed);
            gr_true.second->SetLineStyle(1);
            gr_true.second->SetLineWidth(2);
            cn->cd(1);
            gr_true.first->Draw("same lp");
            cn->cd(2);
            gr_true.second->Draw("same lp");
            legend->AddEntry(gr_true.first, "MC hits");
            // Draw legend
            cn->cd(1);
            legend->Draw();


            // output
            cn->Write(TString::Format("event_%d_track_%d_positions", (int) nevents, trackRow+1).Data());
            // delete graphMap
            for (auto& [k, g] : graphMap) {
                delete g.first;
                delete g.second;
            }
            ///////////////////////////////
            // Display : error variances
            // ////////////////////////////
            TCanvas * cm = new TCanvas(TString::Format("c%d_%d_var", (int) nevents, trackRow+1).Data(),TString::Format("c%d_%d ErrorCovaraince", (int) nevents, trackRow+1).Data(), 1600, 900);
            cm->Divide(3,2);
            std::vector<const char *> name = {"x", "y", "z", "px", "py", "pz"};
            std::vector<const char *> unit = {"mm", "mm", "mm", "MeV", "MeV", "MeV"};
            for (int i = 0; i < 6; i++) {
                cm->cd(i+1);
                graphVariance[0][i]->SetMarkerColor(kRed);
                graphVariance[0][i]->SetMarkerStyle(25);
                graphVariance[0][i]->SetMarkerSize(2);
                graphVariance[1][i]->SetMarkerColor(kBlue);
                graphVariance[1][i]->SetMarkerStyle(26);
                graphVariance[1][i]->SetMarkerSize(2);
                TMultiGraph *mg1 = new TMultiGraph();
                mg1->SetTitle(TString::Format("Variance %s in %s^{2}; Niter; Var(%s)", name[i], unit[i], name[i]).Data());
                mg1->Add(graphVariance[2][i], "l");
                mg1->Add(graphVariance[1][i], "p");
                mg1->Add(graphVariance[0][i], "p");
                mg1->Draw("a");
                TLegend* legend2 = new TLegend(0.1,0.7,0.35,0.9);
                //legend2->AddEntry(graphVariance[2][i], "all");
                legend2->AddEntry(graphVariance[0][i], "prediction", "p");
                legend2->AddEntry(graphVariance[1][i], "correction", "p");
                legend2->Draw();
            }
            // write output
            cm->Write(TString::Format("event_%d_track_%d_variances", (int) nevents, trackRow+1).Data());
            ///////////////////////////////
            // Display : error values
            // ////////////////////////////
            TCanvas * co = new TCanvas(TString::Format("c%d_%d_val", (int) nevents, trackRow+1).Data(),TString::Format("c%d_%d State or Value", (int) nevents, trackRow+1).Data(), 1600, 900);
            co->Divide(3,2);
            for (int i = 0; i < 6; i++) {
                co->cd(i+1);
                graphValue[0][i]->SetMarkerColor(kRed);
                graphValue[0][i]->SetMarkerStyle(25);
                graphValue[0][i]->SetMarkerSize(2);
                graphValue[1][i]->SetMarkerColor(kBlue);
                graphValue[1][i]->SetMarkerStyle(26);
                graphValue[1][i]->SetMarkerSize(2);
                graphValue[3][i]->SetLineColor(kGreen);
                graphValue[3][i]->SetLineWidth(2);
                TMultiGraph *mg1 = new TMultiGraph();
                mg1->SetTitle(TString::Format("Value %s in %s; Niter; %s", name[i], unit[i], name[i]).Data());
                mg1->Add(graphValue[3][i], "l");
                mg1->Add(graphValue[2][i], "l");
                mg1->Add(graphValue[1][i], "p");
                mg1->Add(graphValue[0][i], "p");
                mg1->Draw("a");
                TLegend* legend3 = new TLegend(0.1,0.7,0.35,0.9);
                legend3->AddEntry(graphValue[3][i], "true");
                //legend3->AddEntry(graphValue[2][i], "all");
                legend3->AddEntry(graphValue[0][i], "prediction", "p");
                legend3->AddEntry(graphValue[1][i], "correction", "p");
                legend3->Draw();
            }
            // write output
            co->Write(TString::Format("event_%d_track_%d_values", (int) nevents, trackRow+1).Data());
        }
    }
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output

    f->Close();
    printf("File created : %s\n", output);

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
