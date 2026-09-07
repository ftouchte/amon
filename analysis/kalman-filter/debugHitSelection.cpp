/***********************************************
 * Debug hit selection for 10.6 GeV run
 *
 * @author Felix Touchte Codjo
 * @date September 07, 2026
 * ********************************************/

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

#include "fOptions.h"
#include "debugHitSelection.h"
#include "AhdcCCDB.h"


int main(int argc, char const *argv[]) {

    // record start time
    auto start = std::chrono::high_resolution_clock::now();

    fOptions OPT({"-i", "-o"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::string filename = OPT.GetValue("-i");
    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    if (filename.compare("") == 0 || output.compare("") == 0) {
        printf("Please provide options... (if you specify -i, you should also specify -o)\n");
        return 1;
    }

    // Histograms
    Histograms* histos = new Histograms();

    // Nb events
    long unsigned int nevents = 0;

    int nfile = 0;
    for (auto file : filenames) {
        nfile++;
        printf("> Open file %d/%d : %s\n", nfile, (int) filenames.size(), file.c_str());
        hipo::reader  reader(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        // bank definition
        hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
        hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
        hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
        hipo::bank  trackBank(factory.getSchema("AHDC::track"));
        hipo::bank  runBank(factory.getSchema("RUN::config"));
        hipo::bank  recEventBank(factory.getSchema("REC::Event"));
        hipo::event event;
        long unsigned int nevents_per_file =0;

        // CCDB connection
        AhdcCCDB* ahdcConstants = new AhdcCCDB();
        
        /////////////////////////
        // Loop over events
        /////////////////////////
        while( reader.next()){
            nevents++;
            nevents_per_file++;

            // display progress Bar
            if ((nevents_per_file % 1000 == 0) || ((int) nevents_per_file == reader.getEntries())) {
                progressBar(100.0*nevents_per_file/reader.getEntries());
            }

            // load bank content for this event
            reader.read(event);
            event.getStructure(adcBank);
            event.getStructure(wfBank);
            event.getStructure(hitBank);
            event.getStructure(trackBank);
            event.getStructure(runBank);
            event.getStructure(recEventBank);

            if (nevents == 1) {
                int run = runBank.getInt("run", 0);
                ahdcConstants->setRunNumber(run);
                ahdcConstants->loadConstants();
            }
            

            double startTime = recEventBank.getFloat("startTime", 0);
            histos->H1_startTime->Fill(startTime);

            // any hits
            for (int i = 0; i < adcBank.getRows(); i++) {

                int layer = adcBank.get("layer", i);
                int component = adcBank.get("component", i);

                double adc            = adcBank.getInt("ADC", i);
                double leadingEdgeTime   = adcBank.getFloat("leadingEdgeTime", i);
                double tot = adcBank.getFloat("timeOverThreshold", i);
                double ped         = adcBank.getFloat("ped", i);
                int    wfType            = adcBank.getShort("wfType", i);
                double t0         = ahdcConstants->get_t0(1, layer, component).t0;
                double time = leadingEdgeTime - t0 - startTime;
                
                histos->H1_any_hit_occupancy->Fill(slc2wire(1,layer,component));
                histos->H1_any_hit_amplitude->Fill(adc);
                histos->H1_any_hit_tot->Fill(tot);
                histos->H1_any_hit_time->Fill(time);
                histos->H1_any_hit_ped->Fill(ped);
                histos->H1_any_hit_wfType->Fill(wfType);
            }

            // selected hits
            // not really relevant as we only save hit associated to tracks
            // put here for principle
            // it is also a bug in the reconstruction
            for (int h = 0; h < hitBank.getRows(); h++) {

                int i = hitBank.getShort("id", h);
                
                int layer = adcBank.get("layer", i);
                int component = adcBank.get("component", i);

                double adc            = adcBank.getInt("ADC", i);
                double leadingEdgeTime   = adcBank.getFloat("leadingEdgeTime", i);
                double tot = adcBank.getFloat("timeOverThreshold", i);
                double ped         = adcBank.getFloat("ped", i);
                int    wfType            = adcBank.getShort("wfType", i);
                double t0         = ahdcConstants->get_t0(1, layer, component).t0;
                double time = leadingEdgeTime - t0 - startTime;
                
                histos->H1_selected_hit_occupancy->Fill(slc2wire(1,layer,component));
                histos->H1_selected_hit_amplitude->Fill(adc);
                histos->H1_selected_hit_tot->Fill(tot);
                histos->H1_selected_hit_time->Fill(time);
                histos->H1_selected_hit_ped->Fill(ped);
                histos->H1_selected_hit_wfType->Fill(wfType);

            }

            for (int t = 0; t < trackBank.getRows(); t++) {
                int trackid = trackBank.getInt("trackid", t);
                
                int nhits = trackBank.getInt("n_hits", t);
                histos->H1_track_nhits->Fill(nhits);

                // selected hits
                for (int h = 0; h < hitBank.getRows(); h++) {

                    if (hitBank.getInt("trackid", h) != trackid) continue;

                    int i = hitBank.getShort("id", h) - 1;
                    
                    int layer = adcBank.get("layer", i);
                    int component = adcBank.get("component", i);

                    double adc            = adcBank.getInt("ADC", i);
                    double leadingEdgeTime   = adcBank.getFloat("leadingEdgeTime", i);
                    double tot = adcBank.getFloat("timeOverThreshold", i);
                    double ped         = adcBank.getFloat("ped", i);
                    int    wfType            = adcBank.getShort("wfType", i);
                    double t0         = ahdcConstants->get_t0(1, layer, component).t0;
                    double time = leadingEdgeTime - t0 - startTime;
                    
                    histos->H1_track_hit_occupancy->Fill(slc2wire(1,layer,component));
                    histos->H1_track_hit_amplitude->Fill(adc);
                    histos->H1_track_hit_tot->Fill(tot);
                    histos->H1_track_hit_time->Fill(time);
                    histos->H1_track_hit_ped->Fill(ped);
                    histos->H1_track_hit_wfType->Fill(wfType);
                    

                }
            }
            
        } // loop over events 

    } // loop over files

    // occupancy renormalization
    histos->H1_any_hit_occupancy->Scale(100.0/nevents);
    histos->H1_selected_hit_occupancy->Scale(100.0/nevents);
    histos->H1_track_hit_occupancy->Scale(100.0/nevents);

    TFile *f = new TFile(output.c_str(), "RECREATE");

    histos->WriteIn(f);

    f->Close();


    // end of the program
    printf("* nevents : %ld\n", nevents);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("* time elapsed : %lf seconds\n", elapsed.count());
    return 0;

}


/**
 * @brief Convert (sector, layer, component) to a unqiue wire id (number betwwen 0 and 575)
 * 
 * @param sector (not used)
 * @param layer 
 * @param component 
 * @return unique wire id
 */
int slc2wire(int sector, int layer, int component) {
    if (layer == 11) {
        return component - 1;
    } 
    else if (layer == 21) {
        return 47 + component - 1;
    } 
    else if (layer == 22) {
        return 47 + 56 + component - 1;
    } 
    else if (layer == 31) {
        return 47 + 56 + 56 + component - 1;
    } 
    else if (layer == 32) {
        return 47 + 56 + 56 + 72 + component - 1;
    } 
    else if (layer == 41) {
        return 47 + 56 + 56 + 72 + 72 + component - 1;
    } 
    else if (layer == 42) {
        return 47 + 56 + 56 + 72 + 72 + 87 + component - 1;
    } 
    else if (layer == 51) {
        return 47 + 56 + 56 + 72 + 72 + 87 + 87 + component - 1;
    } else {
        return -1; // not a ahdc wire
    }
}

int layer2number(int digit) {
    if      (digit == 11) {
        return 1;
    } 
    else if (digit == 21) {
        return 2;
    } 
    else if (digit == 22) {
        return 3;
    } 
    else if (digit == 31) {
        return 4;
    } 
    else if (digit == 32) {
        return 5;
    } 
    else if (digit == 41) {
        return 6;
    } 
    else if (digit == 42) {
        return 7;
    } 
    else if (digit == 51) {
        return 8;
    } else {
        return 0; // not a layer, can encode all layers
    }
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