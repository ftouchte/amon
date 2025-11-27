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
void computeSphericalVariance(double mu_x, double mu_y, double mu_z, double var_x, double var_y, double var_z, double & var_r, double & var_theta, double & var_phi, const char * nature = "");


int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
   const char * output = "./output/kfmon-v31.root";
    
    /////////////////////////
    /// simu
    /// /////////////////////
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v31.hipo";
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
    hipo::bank  runConfigBank(factory.getSchema("RUN::config"));
    hipo::bank  recBank(factory.getSchema("REC::Particle"));
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

    // Histograms
    TH2D* H2_corr_phi = new TH2D("corr_phi", "#phi_{mc} vs #phi_{track}; #phi_{track} (deg); #phi_{mc} (deg)", 100, 0, 361, 100, 0, 361); 
    TH2D* H2_corr_theta = new TH2D("corr_theta", "#theta_{mc} vs #theta_{track}; #theta_{track} (deg); #theta_{mc} (deg)", 100, 0, 181, 100, 50, 130); 
    TH2D* H2_corr_p = new TH2D("corr_p", "p_{mc} vs p_{track}; p_{track} (MeV); p_{mc} (MeV)", 100, 190, 310, 100, 190, 310); 
    TH2D* H2_corr_vz = new TH2D("corr_vz", "vz_{mc} vs vz_{track}; vz_{track} (cm); vz_{mc} (cm)", 100, -16, 16, 100, -16, 16); 
    TH1D* H1_delta_phi = new TH1D("delta_phi", "#Delta #phi = #phi_{mc} - #phi_{track}; #Delta #phi (deg); #count", 100, -12, 12);
    TH1D* H1_delta_vz = new TH1D("delta_vz", "#Delta vz = vz_{mc} - vz_{track}; #Delta vz (cm); #count", 100, -10, 10);
    TH1D* H1_residual = new TH1D("residual", "residual; residual (mm); #count", 100, -2, 2);
    TH1D* H1_residual_filtered = new TH1D("residual_filtered", "residual; residual (mm); #count", 100, -2, 2);
    TH1D* H1_time = new TH1D("time", "time; time (ns); #count", 100, 0, 250);
    TH1D* H1_distance = new TH1D("distance", "distance; distance (mm); #count", 100, 0, 4);
    TH1I* H1_wfType = new TH1I("wfType", "wfType; wfType; #count", 7, 0, 7);

    // Loop over events
    while( reader.next()){

        //if (nevents > 10) break;
        nevents++;
        // Progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        event.getStructure(trackBank);
        event.getStructure(kfmonBank);
        event.getStructure(track0Bank);
        event.getStructure(hitBank);
        event.getStructure(mcBank);
        event.getStructure(trueBank);
        event.getStructure(runConfigBank);
        event.getStructure(recBank);
        if (recBank.getRows() < 1) continue; 


        nMCtracks += mcBank.getRows();
        nKFtracks += trackBank.getRows();
        { // Histograms
            if (trackBank.getRows() < 1) continue;
            // mc
            double mc_px = mcBank.get("px", 0)*1000; // GeV to MeV
            double mc_py = mcBank.get("py", 0)*1000;
            double mc_pz = mcBank.get("pz", 0)*1000;
            double mc_p;
            double mc_theta;
            double mc_phi;
            futils::cart2polarDEG(mc_px,mc_py,mc_pz,mc_p,mc_theta,mc_phi);
            // track
            double track_px = trackBank.get("px", 0); // MeV
            double track_py = trackBank.get("py", 0);
            double track_pz = trackBank.get("pz", 0);
            double track_p;
            double track_theta;
            double track_phi;
            futils::cart2polarDEG(track_px,track_py,track_pz,track_p,track_theta,track_phi);
            double mc_vz = mcBank.get("vz", 0); // cm
            double track_vz = trackBank.get("z", 0)*0.1; // mm to cm
            // Fill historgrams
            H2_corr_phi->Fill(track_phi, mc_phi);
            H2_corr_theta->Fill(track_theta, mc_theta);
            H2_corr_p->Fill(track_p, mc_p);
            H2_corr_vz->Fill(track_vz, mc_vz);
            H1_delta_phi->Fill(mc_phi-track_phi);
            H1_delta_vz->Fill(mc_vz-track_vz);
            for (int i = 0; i < hitBank.getRows(); i++) {
                if ((int) hitBank.get("trackid", i) == (int) trackBank.get("trackid", 0)) {
                    H1_residual->Fill(hitBank.get("residual", i));
                    if (std::abs(hitBank.get("residual", i)) > 1e-9) {
                        H1_residual_filtered->Fill(hitBank.get("residual", i));
                    }
                    H1_time->Fill(hitBank.get("time", i));
                    H1_distance->Fill(hitBank.get("doca", i));
                    H1_wfType->Fill(adcBank.get("wfType", hitBank.get("id", i) - 1));
                }
            }
        }
        
        if (nevents > 10) continue;
        int evt_number = runConfigBank.get("event", 0);
        //if (evt_number < 4501 || evt_number > 4510) continue;
        //printf("evt n° : %d\n", (int) runConfigBank.get("event", 0));
        // loop over track
        for (int trackRow = 0; trackRow < trackBank.getRows(); trackRow++) {
            ///////////////////////////////
            // Display : states/positions
            // ////////////////////////////
            TCanvas * cn = new TCanvas(TString::Format("c%d_%d", (int) evt_number, trackRow+1).Data(),TString::Format("c%d_%d ahdcView", (int) evt_number, trackRow+1).Data(), 1600, 900);
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
            // Variances of x, y, z, modulus, theta, phi (spherical coordinates) around the beam axis (indicator ==0) 
            //std::vector<TGraph*> graphVariance = {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()};
            // Position
            std::vector<std::vector<TGraph*>> graphPosVariance = {
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // prediction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // correction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}  // both of them
            };
            // Momentum
            std::vector<std::vector<TGraph*>> graphMomVariance = {
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // prediction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // correction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}  // both of them
            };
            std::vector<std::vector<TGraph*>> graphPosValue = {
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // prediction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}, // correction
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()},  // both of them
                {new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph(), new TGraph()}  // mc value
            };
            std::vector<std::vector<TGraph*>> graphMomValue = {
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
                double var_pos_r;
                double var_pos_theta;
                double var_pos_phi;
                double var_mom_p;
                double var_mom_theta;
                double var_mom_phi;
                computeSphericalVariance(x, y, z, var_x, var_y, var_z, var_pos_r, var_pos_theta, var_pos_phi, "position");
                computeSphericalVariance(px, py, pz, var_px, var_py, var_pz, var_mom_p, var_mom_theta, var_mom_phi, "momentum");
                double pos_r;
                double pos_theta;
                double pos_phi;
                double mom_p;
                double mom_theta;
                double mom_phi;
                futils::cart2polarDEG(x,y,z,pos_r,pos_theta,pos_phi);
                futils::cart2polarDEG(px,py,pz,mom_p,mom_theta,mom_phi);
                double mc_vx = mcBank.get("vx", 0)*10; // cm to mm 
                double mc_vy = mcBank.get("vy", 0)*10; 
                double mc_vz = mcBank.get("vz", 0)*10; 
                double mc_px = mcBank.get("px", 0)*1000; // GeV to MeV
                double mc_py = mcBank.get("py", 0)*1000;
                double mc_pz = mcBank.get("pz", 0)*1000;
                double mc_pos_r;
                double mc_pos_theta;
                double mc_pos_phi;
                double mc_mom_p;
                double mc_mom_theta;
                double mc_mom_phi;
                futils::cart2polarDEG(mc_vx,mc_vy,mc_vz,mc_pos_r,mc_pos_theta,mc_pos_phi);
                futils::cart2polarDEG(mc_px,mc_py,mc_pz,mc_mom_p,mc_mom_theta,mc_mom_phi);

                // it : iterator, is an object that points to a std::pair containing 
                // it->first is the key, here the std::tuple
                // it->second is the value, here the TGraph*
                // inserted is a boolean, is false if the element already exist and true if it is a new element
                // if false, "it" corresponds to the existing element
                ///////////////////////////////// 
                /// cut on Niter
                //if (Niter >= 10) continue;
                /////////////////////////////////
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
                       graphPosVariance[2][0]->AddPoint(Niter, var_x); 
                       graphPosVariance[2][1]->AddPoint(Niter, var_y); 
                       graphPosVariance[2][2]->AddPoint(Niter, var_z); 
                       graphPosVariance[2][3]->AddPoint(Niter, var_pos_r); 
                       graphPosVariance[2][4]->AddPoint(Niter, var_pos_theta); 
                       graphPosVariance[2][5]->AddPoint(Niter, var_pos_phi); 
                       graphMomVariance[2][0]->AddPoint(Niter, var_px); 
                       graphMomVariance[2][1]->AddPoint(Niter, var_py); 
                       graphMomVariance[2][2]->AddPoint(Niter, var_pz);
                       graphMomVariance[2][3]->AddPoint(Niter, var_mom_p); 
                       graphMomVariance[2][4]->AddPoint(Niter, var_mom_theta); 
                       graphMomVariance[2][5]->AddPoint(Niter, var_mom_phi);
                       if (status == 0) { // prediction
                           graphPosVariance[0][0]->AddPoint(Niter, var_x); 
                           graphPosVariance[0][1]->AddPoint(Niter, var_y); 
                           graphPosVariance[0][2]->AddPoint(Niter, var_z); 
                           graphPosVariance[0][3]->AddPoint(Niter, var_pos_r); 
                           graphPosVariance[0][4]->AddPoint(Niter, var_pos_theta); 
                           graphPosVariance[0][5]->AddPoint(Niter, var_pos_phi); 
                           graphMomVariance[0][0]->AddPoint(Niter, var_px); 
                           graphMomVariance[0][1]->AddPoint(Niter, var_py); 
                           graphMomVariance[0][2]->AddPoint(Niter, var_pz);
                           graphMomVariance[0][3]->AddPoint(Niter, var_mom_p); 
                           graphMomVariance[0][4]->AddPoint(Niter, var_mom_theta); 
                           graphMomVariance[0][5]->AddPoint(Niter, var_mom_phi);
                       } 
                       else { // correction, status == 1
                           graphPosVariance[1][0]->AddPoint(Niter, var_x); 
                           graphPosVariance[1][1]->AddPoint(Niter, var_y); 
                           graphPosVariance[1][2]->AddPoint(Niter, var_z); 
                           graphPosVariance[1][3]->AddPoint(Niter, var_pos_r); 
                           graphPosVariance[1][4]->AddPoint(Niter, var_pos_theta); 
                           graphPosVariance[1][5]->AddPoint(Niter, var_pos_phi); 
                           graphMomVariance[1][0]->AddPoint(Niter, var_px); 
                           graphMomVariance[1][1]->AddPoint(Niter, var_py); 
                           graphMomVariance[1][2]->AddPoint(Niter, var_pz);
                           graphMomVariance[1][3]->AddPoint(Niter, var_mom_p); 
                           graphMomVariance[1][4]->AddPoint(Niter, var_mom_theta); 
                           graphMomVariance[1][5]->AddPoint(Niter, var_mom_phi);
                       }
                    }
                }
                { // Value graph
                    if (indicator == 0) { // beam level
                       graphPosValue[2][0]->AddPoint(Niter, x); 
                       graphPosValue[2][1]->AddPoint(Niter, y); 
                       graphPosValue[2][2]->AddPoint(Niter, z); 
                       graphPosValue[2][3]->AddPoint(Niter, pos_r); 
                       graphPosValue[2][4]->AddPoint(Niter, pos_theta); 
                       graphPosValue[2][5]->AddPoint(Niter, pos_phi);
                       graphMomValue[2][0]->AddPoint(Niter, px); 
                       graphMomValue[2][1]->AddPoint(Niter, py); 
                       graphMomValue[2][2]->AddPoint(Niter, pz); 
                       graphMomValue[2][3]->AddPoint(Niter, mom_p); 
                       graphMomValue[2][4]->AddPoint(Niter, mom_theta); 
                       graphMomValue[2][5]->AddPoint(Niter, mom_phi);
                       if (status == 0) { // prediction
                           graphPosValue[0][0]->AddPoint(Niter, x); 
                           graphPosValue[0][1]->AddPoint(Niter, y); 
                           graphPosValue[0][2]->AddPoint(Niter, z); 
                           graphPosValue[0][3]->AddPoint(Niter, pos_r); 
                           graphPosValue[0][4]->AddPoint(Niter, pos_theta); 
                           graphPosValue[0][5]->AddPoint(Niter, pos_phi); 
                           graphMomValue[0][0]->AddPoint(Niter, px); 
                           graphMomValue[0][1]->AddPoint(Niter, py); 
                           graphMomValue[0][2]->AddPoint(Niter, pz);
                           graphMomValue[0][3]->AddPoint(Niter, mom_p); 
                           graphMomValue[0][4]->AddPoint(Niter, mom_theta); 
                           graphMomValue[0][5]->AddPoint(Niter, mom_phi);
                       } 
                       else { // correction, status == 1
                           graphPosValue[1][0]->AddPoint(Niter, x); 
                           graphPosValue[1][1]->AddPoint(Niter, y); 
                           graphPosValue[1][2]->AddPoint(Niter, z); 
                           graphPosValue[1][3]->AddPoint(Niter, pos_r); 
                           graphPosValue[1][4]->AddPoint(Niter, pos_theta); 
                           graphPosValue[1][5]->AddPoint(Niter, pos_phi); 
                           graphMomValue[1][0]->AddPoint(Niter, px); 
                           graphMomValue[1][1]->AddPoint(Niter, py); 
                           graphMomValue[1][2]->AddPoint(Niter, pz);
                           graphMomValue[1][3]->AddPoint(Niter, mom_p); 
                           graphMomValue[1][4]->AddPoint(Niter, mom_theta); 
                           graphMomValue[1][5]->AddPoint(Niter, mom_phi);
                       }
                       graphPosValue[3][0]->AddPoint(Niter, mc_vx); 
                       graphPosValue[3][1]->AddPoint(Niter, mc_vy); 
                       graphPosValue[3][2]->AddPoint(Niter, mc_vz); 
                       graphPosValue[3][3]->AddPoint(Niter, mc_pos_r); 
                       graphPosValue[3][4]->AddPoint(Niter, mc_pos_theta); 
                       graphPosValue[3][5]->AddPoint(Niter, mc_pos_phi); 
                       graphMomValue[3][0]->AddPoint(Niter, mc_px); 
                       graphMomValue[3][1]->AddPoint(Niter, mc_py); 
                       graphMomValue[3][2]->AddPoint(Niter, mc_pz);
                       graphMomValue[3][3]->AddPoint(Niter, mc_mom_p); 
                       graphMomValue[3][4]->AddPoint(Niter, mc_mom_theta); 
                       graphMomValue[3][5]->AddPoint(Niter, mc_mom_phi);
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
            cn->Write(TString::Format("event_%d_track_%d_positions", (int) evt_number, trackRow+1).Data());
            // delete graphMap
            for (auto& [k, g] : graphMap) {
                delete g.first;
                delete g.second;
            }
            ///////////////////////////////
            // Display : error variances
            // ////////////////////////////
            TCanvas * cvar = new TCanvas(TString::Format("c%d_%d_var_pos", (int) evt_number, trackRow+1).Data(),TString::Format("c%d_%d ErrorVariance", (int) evt_number, trackRow+1).Data(), 1600, 900);
            cvar->Divide(4,3);
            std::vector<const char *> name_pos = {"Position x", "Position y", "Position z", "Position R", "Position theta", "Position phi"};
            std::vector<const char *> unit_pos = {"mm", "mm", "mm", "mm", "deg", "deg"};
            std::vector<const char *> name_mom = {"Momentum px", "Momentum py", "Momentum pz", "Momentum p", "Momentum theta", "Momentum phi"};
            std::vector<const char *> unit_mom = {"MeV", "MeV", "MeV", "MeV", "deg", "deg"};
            for (int i = 0; i < 6; i++) {
                {// Pos
                    cvar->cd(i+1);
                    graphPosVariance[0][i]->SetMarkerColor(kRed);
                    graphPosVariance[0][i]->SetMarkerStyle(25);
                    graphPosVariance[0][i]->SetMarkerSize(2);
                    graphPosVariance[1][i]->SetMarkerColor(kBlue);
                    graphPosVariance[1][i]->SetMarkerStyle(26);
                    graphPosVariance[1][i]->SetMarkerSize(2);
                    TMultiGraph *mg1 = new TMultiGraph();
                    mg1->SetTitle(TString::Format("Variance %s in %s^{2}; Niter; Var(%s)", name_pos[i], unit_pos[i], name_pos[i]).Data());
                    mg1->Add(graphPosVariance[2][i], "l");
                    mg1->Add(graphPosVariance[1][i], "p");
                    mg1->Add(graphPosVariance[0][i], "p");
                    TLegend* legend1 = new TLegend(0.1,0.7,0.35,0.9);
                    legend1->AddEntry(graphPosVariance[0][i], "prediction", "p");
                    legend1->AddEntry(graphPosVariance[1][i], "correction", "p");
                    mg1->Draw("a");
                    legend1->Draw();
                }
                {// Mom
                    cvar->cd(6+i+1);
                    graphMomVariance[0][i]->SetMarkerColor(kRed);
                    graphMomVariance[0][i]->SetMarkerStyle(25);
                    graphMomVariance[0][i]->SetMarkerSize(2);
                    graphMomVariance[1][i]->SetMarkerColor(kBlue);
                    graphMomVariance[1][i]->SetMarkerStyle(26);
                    graphMomVariance[1][i]->SetMarkerSize(2);
                    TMultiGraph *mg2 = new TMultiGraph();
                    mg2->SetTitle(TString::Format("Variance %s in %s^{2}; Niter; Var(%s)", name_mom[i], unit_mom[i], name_mom[i]).Data());
                    mg2->Add(graphMomVariance[2][i], "l");
                    mg2->Add(graphMomVariance[1][i], "p");
                    mg2->Add(graphMomVariance[0][i], "p");
                    TLegend* legend2 = new TLegend(0.1,0.7,0.35,0.9);
                    legend2->AddEntry(graphMomVariance[0][i], "prediction", "p");
                    legend2->AddEntry(graphMomVariance[1][i], "correction", "p");
                    mg2->Draw("a");
                    legend2->Draw();
                }
            }
            // write output
            cvar->Write(TString::Format("event_%d_track_%d_variances", (int) evt_number, trackRow+1).Data());
            ///////////////////////////////
            // Display : error values
            // ////////////////////////////
            TCanvas * cval = new TCanvas(TString::Format("c%d_%d_val", (int) evt_number, trackRow+1).Data(),TString::Format("c%d_%d Estimation", (int) evt_number, trackRow+1).Data(), 1600, 900);
            cval->Divide(4,3);
            for (int i = 0; i < 6; i++) {
                { // Pos
                    cval->cd(i+1);
                    graphPosValue[0][i]->SetMarkerColor(kRed);
                    graphPosValue[0][i]->SetMarkerStyle(25);
                    graphPosValue[0][i]->SetMarkerSize(2);
                    graphPosValue[1][i]->SetMarkerColor(kBlue);
                    graphPosValue[1][i]->SetMarkerStyle(26);
                    graphPosValue[1][i]->SetMarkerSize(2);
                    graphPosValue[3][i]->SetLineColor(kGreen);
                    graphPosValue[3][i]->SetLineWidth(2);
                    TMultiGraph *mg1 = new TMultiGraph();
                    mg1->SetTitle(TString::Format("Estimation %s in %s; Niter; %s", name_pos[i], unit_pos[i], name_pos[i]).Data());
                    mg1->Add(graphPosValue[3][i], "l");
                    mg1->Add(graphPosValue[2][i], "l");
                    mg1->Add(graphPosValue[1][i], "p");
                    mg1->Add(graphPosValue[0][i], "p");
                    mg1->Draw("a");
                    TLegend* legend3 = new TLegend(0.1,0.7,0.35,0.9);
                    legend3->AddEntry(graphPosValue[3][i], "true");
                    legend3->AddEntry(graphPosValue[0][i], "prediction", "p");
                    legend3->AddEntry(graphPosValue[1][i], "correction", "p");
                    legend3->Draw();
                }
                { // Mom
                    cval->cd(6+i+1);
                    graphMomValue[0][i]->SetMarkerColor(kRed);
                    graphMomValue[0][i]->SetMarkerStyle(25);
                    graphMomValue[0][i]->SetMarkerSize(2);
                    graphMomValue[1][i]->SetMarkerColor(kBlue);
                    graphMomValue[1][i]->SetMarkerStyle(26);
                    graphMomValue[1][i]->SetMarkerSize(2);
                    graphMomValue[3][i]->SetLineColor(kGreen);
                    graphMomValue[3][i]->SetLineWidth(2);
                    TMultiGraph *mg1 = new TMultiGraph();
                    mg1->SetTitle(TString::Format("Estimation %s in %s; Niter; %s", name_mom[i], unit_mom[i], name_mom[i]).Data());
                    mg1->Add(graphMomValue[3][i], "l");
                    mg1->Add(graphMomValue[2][i], "l");
                    mg1->Add(graphMomValue[1][i], "p");
                    mg1->Add(graphMomValue[0][i], "p");
                    mg1->Draw("a");
                    TLegend* legend3 = new TLegend(0.1,0.7,0.35,0.9);
                    legend3->AddEntry(graphMomValue[3][i], "true");
                    legend3->AddEntry(graphMomValue[0][i], "prediction", "p");
                    legend3->AddEntry(graphMomValue[1][i], "correction", "p");
                    legend3->Draw();
                }
            }
            // write output
            cval->Write(TString::Format("event_%d_track_%d_values", (int) evt_number, trackRow+1).Data());
        }
    }
    printf("nevents    : %ld \n", nevents);
    printf("nMCtrack    : %ld \n", nMCtracks);
    printf("nKFtrack    : %ld \n", nKFtracks);
    printf("percentage  : %.2lf %%\n", 100.0*nKFtracks/nMCtracks);
    // output
    // Write histograms
    H2_corr_phi->Write("corr_phi");
    H2_corr_theta->Write("corr_theta");
    H2_corr_p->Write("corr_p");
    H2_corr_vz->Write("corr_vz");
    H1_delta_phi->Write("delta_phi");
    H1_delta_vz->Write("delta_vz");
    H1_residual->Write("residual");
    H1_residual_filtered->Write("residual_filtered");
    H1_time->Write("time");
    H1_distance->Write("distance");
    H1_wfType->Write("wfType");


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

void computeSphericalVariance(double mu_x, double mu_y, double mu_z, double var_x, double var_y, double var_z, double & var_r, double & var_theta, double & var_phi, const char * nature) {
    // prevent NaN ot Inf values
    // if mu_x and mu_y are strictly equal to 0 (at beginning of the Kalman Filter for example)
    // phi can not be defined
    // theta is equals to 0 but its derivative is infinite at 0
    // we then assume that all the variance if given by the z component to the r component
    // I assume x, y and z independents but in reality we can correlations
    if ((std::abs(mu_x) < 1e-9) && (std::abs(mu_y) < 1e-9)) {
        var_r = var_z;
        var_theta = 0;
        var_phi = 0;
    } else {
        // var_r
        double r = sqrt(mu_x*mu_x+mu_y*mu_y+mu_z*mu_z);
        double drdx = mu_x/r;
        double drdy = mu_y/r;
        double drdz = mu_z/r;
        var_r = drdx*drdx*var_x + drdy*drdy*var_y + drdz*drdz*var_z;
        // var_theta
        double dthetadz = (-1/r)/sqrt(1-pow(mu_z/r,2.0));
        double dthetady = (1/r)*(mu_z/r)*(mu_y/r)/sqrt(1-pow(mu_z/r,2.0));
        double dthetadx = (1/r)*(mu_z/r)*(mu_x/r)/sqrt(1-pow(mu_z/r,2.0));
        var_theta = dthetadx*dthetadx*var_x + dthetady*dthetady*var_y + dthetadz*dthetadz*var_z;
        var_theta *= pow(180.0/M_PI,2.0); // convert to deg²
        // var_phi
        double dphidz = 0;
        double rho2 = mu_x*mu_x+mu_y*mu_y;
        double dphidy = mu_x/rho2;
        double dphidx = -mu_y/rho2;
        var_phi = dphidx*dphidx*var_x + dphidy*dphidy*var_y + dphidz*dphidz*var_z;
        var_phi *= pow(180.0/M_PI,2.0); // convert to deg²
    }
}

