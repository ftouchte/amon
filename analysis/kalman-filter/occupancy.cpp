/***********************************************
 * Occupancy study
 * 
 * Study the effect of the wfType cut vs the ADC cut.
 * 
 * To do so, one only need to run the HitReader inside the AHDCEngine.
 * 
 * In the current state of the reconstruction, at the end of the ALERTEngine, only hits associated with a tracks are recorded. That why I only need the HitReader.
 *
 * @author Felix Touchte Codjo
 * @date September 02, 2026
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


int main(int argc, char const *argv[]) {

    // record start time
    auto start = std::chrono::high_resolution_clock::now();

    fOptions OPT({"-i", "-o", "-v", "-simu"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::string filename = OPT.GetValue("-i");
    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");
    std::string version = OPT.GetValue("-v");
    bool IsMC = OPT.GetValue("-simu").compare("true") == 0;

    int version_number = 1'000'000;

    if (version.compare("") != 0) {
        filename = std::string("/home/touchte-codjo/Desktop/hipofiles/kalman-filter/rec-data-r22712-v") + version + ".hipo"; 
        output = std::string("./output/kfmon_data_r22712_v") + version + ".root";
        version_number = std::atoi(version.c_str());
        if (IsMC) {
            filename = std::string("/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v") + version + ".hipo"; 
            output = std::string("./output/kfmon_data_rsimu_v") + version + ".root";
        }
        filenames = {filename};
    } else {
        if (filename.compare("") == 0 || output.compare("") == 0) {
            printf("Please provide options... (if you specify -i, you should also specify -o)\n");
            return 1;
        }
    }

    // Histograms


    TH1D* H1_selected_hit_occupancy = new TH1D("occupancy_from_selected_hits", "occupancy; wire; occupancy [%]", 576, 0, 576); // after HitReader
    TH1D* H1_all_hit_occupancy = new TH1D("occupancy_from_all_hits", "occupancy; wire; occupancy [%]", 576, 0, 576); // before HitReader
    TH1D* H1_amplitude = new TH1D("selected_hit_amplitude", "amplitude ; amplitude (ADC); count",  100, 0, 4000);
    TH1I* H1_wfType = new TH1I("selected_hit_wfType", "wfType; wfType; #count", 7, 0, 7);
    TH2D* H2_wire_occupancy = new TH2D("wire_occupancy", "Wire occupancy; wire; layer", 99, 1, 100, 8, 1, 9);


    long unsigned int nevents =0;
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
        hipo::event event;
        long unsigned int nevents_per_file =0;
        
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

            // occupancy
            for (int i = 0; i < adcBank.getRows(); i++) {
                int layer = adcBank.get("layer", i);
                int component = adcBank.get("component", i);
                H1_all_hit_occupancy->Fill(slc2wire(1,layer,component));
                H2_wire_occupancy->Fill(component, layer2number(layer));
            }

            for (int i = 0; i < hitBank.getRows(); i++) {
                int layer = 10*hitBank.get("superlayer", i) + hitBank.get("layer", i);
                int component = hitBank.get("wire", i);
                H1_selected_hit_occupancy->Fill(slc2wire(1,layer,component));

                int adcRow = hitBank.get("id", i)-1;
                double adc = adcBank.getInt("ADC", adcRow);
                double wfType = adcBank.getShort("wfType", adcRow);
                H1_amplitude->Fill(adc);
                H1_wfType->Fill(wfType);
            }
            
        } // loop over events 

    } // loop over files

    TFile *f = new TFile(output.c_str(), "RECREATE");

    H1_selected_hit_occupancy->Scale(100.0/nevents);
    H1_all_hit_occupancy->Scale(100.0/nevents);
    H2_wire_occupancy->Scale(100.0/nevents);

    H1_selected_hit_occupancy->Write(H1_selected_hit_occupancy->GetName());
    H1_all_hit_occupancy->Write(H1_all_hit_occupancy->GetName());
    H1_amplitude->Write(H1_amplitude->GetName());
    H1_wfType->Write(H1_wfType->GetName());
    H2_wire_occupancy->Write(H2_wire_occupancy->GetName());

    f->Close();


    // end of the program
    printf("* nevents : %ld\n", nevents);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("* time elapsed : %lf seconds\n", elapsed.count());
    return 0;

}