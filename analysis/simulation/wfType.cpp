/***********************************************
 * Study wfType
 *
 * @author Felix Touchte Codjo
 * @date September 19, 2025
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

#include "AhdcCCDB.h"
#include "futils.h"

// utilities
void progressBar(int state, int bar_length = 100);
// end utilities

int main(int argc, char const *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    
    const char* output = "./output/wfType_study.root";
    const char* filename = "/home/touchte-codjo/Desktop/hipofiles/wfType/new-classification.hipo";
    printf("> filename : %s\n", filename);
    hipo::reader  reader(filename);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
    hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
    hipo::event event;
    long unsigned int nevents =0;
    // Histograms
    std::vector<TH1D*> VecH1_wfType;
    VecH1_wfType.push_back(new TH1D("wfType_wfType=0", "wfType_wfType=0; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=1", "wfType_wfType=1; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=2", "wfType_wfType=2; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=3", "wfType_wfType=3; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=4", "wfType_wfType=4; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=5", "wfType_wfType=5; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType=6", "wfType_wfType=6; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=0", "wfType_wfType<=0; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=1", "wfType_wfType<=1; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=2", "wfType_wfType<=2; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=3", "wfType_wfType<=3; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=4", "wfType_wfType<=4; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=5", "wfType_wfType<=5; wfType (ns); count", 7, 0, 7));
    VecH1_wfType.push_back(new TH1D("wfType_wfType<=6", "wfType_wfType<=6; wfType (ns); count", 7, 0, 7));
    std::vector<TH1D*> VecH1_time;
    VecH1_time.push_back(new TH1D("time_wfType=0", "time_wfType=0; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=1", "time_wfType=1; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=2", "time_wfType=2; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=3", "time_wfType=3; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=4", "time_wfType=4; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=5", "time_wfType=5; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType=6", "time_wfType=6; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=0", "time_wfType<=0; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=1", "time_wfType<=1; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=2", "time_wfType<=2; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=3", "time_wfType<=3; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=4", "time_wfType<=4; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=5", "time_wfType<=5; time (ns); count", 100, 0, 400));
    VecH1_time.push_back(new TH1D("time_wfType<=6", "time_wfType<=6; time (ns); count", 100, 0, 400));
    std::vector<TH1D*> VecH1_leadingEdgeTime;
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=0", "leadingEdgeTime_wfType=0; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=1", "leadingEdgeTime_wfType=1; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=2", "leadingEdgeTime_wfType=2; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=3", "leadingEdgeTime_wfType=3; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=4", "leadingEdgeTime_wfType=4; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=5", "leadingEdgeTime_wfType=5; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType=6", "leadingEdgeTime_wfType=6; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=0", "leadingEdgeTime_wfType<=0; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=1", "leadingEdgeTime_wfType<=1; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=2", "leadingEdgeTime_wfType<=2; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=3", "leadingEdgeTime_wfType<=3; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=4", "leadingEdgeTime_wfType<=4; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=5", "leadingEdgeTime_wfType<=5; leadingEdgeTime (ns); count", 100, 200, 500));
    VecH1_leadingEdgeTime.push_back(new TH1D("leadingEdgeTime_wfType<=6", "leadingEdgeTime_wfType<=6; leadingEdgeTime (ns); count", 100, 200, 500));
    std::vector<TH1D*> VecH1_timeOverThreshold;
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=0", "timeOverThreshold_wfType=0; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=1", "timeOverThreshold_wfType=1; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=2", "timeOverThreshold_wfType=2; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=3", "timeOverThreshold_wfType=3; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=4", "timeOverThreshold_wfType=4; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=5", "timeOverThreshold_wfType=5; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType=6", "timeOverThreshold_wfType=6; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=0", "timeOverThreshold_wfType<=0; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=1", "timeOverThreshold_wfType<=1; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=2", "timeOverThreshold_wfType<=2; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=3", "timeOverThreshold_wfType<=3; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=4", "timeOverThreshold_wfType<=4; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=5", "timeOverThreshold_wfType<=5; timeOverThreshold (ns); count", 100, 150, 650));
    VecH1_timeOverThreshold.push_back(new TH1D("timeOverThreshold_wfType<=6", "timeOverThreshold_wfType<=6; timeOverThreshold (ns); count", 100, 150, 650));
    std::vector<TH1D*> VecH1_adcMax;
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=0", "adcMax_wfType=0; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=1", "adcMax_wfType=1; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=2", "adcMax_wfType=2; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=3", "adcMax_wfType=3; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=4", "adcMax_wfType=4; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=5", "adcMax_wfType=5; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType=6", "adcMax_wfType=6; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=0", "adcMax_wfType<=0; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=1", "adcMax_wfType<=1; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=2", "adcMax_wfType<=2; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=3", "adcMax_wfType<=3; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=4", "adcMax_wfType<=4; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=5", "adcMax_wfType<=5; adcMax (ns); count", 100, 0, 3700));
    VecH1_adcMax.push_back(new TH1D("adcMax_wfType<=6", "adcMax_wfType<=6; adcMax (ns); count", 100, 0, 3700));
    std::vector<TH1D*> VecH1_adcOffset;
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=0", "adcOffset_wfType=0; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=1", "adcOffset_wfType=1; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=2", "adcOffset_wfType=2; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=3", "adcOffset_wfType=3; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=4", "adcOffset_wfType=4; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=5", "adcOffset_wfType=5; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType=6", "adcOffset_wfType=6; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=0", "adcOffset_wfType<=0; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=1", "adcOffset_wfType<=1; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=2", "adcOffset_wfType<=2; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=3", "adcOffset_wfType<=3; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=4", "adcOffset_wfType<=4; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=5", "adcOffset_wfType<=5; adcOffset (ns); count", 100, 100, 400));
    VecH1_adcOffset.push_back(new TH1D("adcOffset_wfType<=6", "adcOffset_wfType<=6; adcOffset (ns); count", 100, 100, 400));
    std::vector<TH1D*> VecH1_occupancy;
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=0", "occupancy_occupancy=0; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=1", "occupancy_occupancy=1; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=2", "occupancy_occupancy=2; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=3", "occupancy_occupancy=3; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=4", "occupancy_occupancy=4; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=5", "occupancy_occupancy=5; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy=6", "occupancy_occupancy=6; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=0", "occupancy_occupancy<=0; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=1", "occupancy_occupancy<=1; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=2", "occupancy_occupancy<=2; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=3", "occupancy_occupancy<=3; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=4", "occupancy_occupancy<=4; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=5", "occupancy_occupancy<=5; occupancy (ns); count [%]", 576, 0, 576));
    VecH1_occupancy.push_back(new TH1D("occupancy_occupancy<=6", "occupancy_occupancy<=6; occupancy (ns); count [%]", 576, 0, 576));

    AhdcCCDB ahdcConstants;
    // Loop over events
    while( reader.next()){
        nevents++;
        // Progress Bar
        if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
            progressBar(100.0*nevents/reader.getEntries());
        }
        reader.read(event);
        event.getStructure(adcBank);
        event.getStructure(wfBank);
        
        // Loop over waveforms
        for (int i = 0; i < adcBank.getRows(); i++) {
            double sector    = adcBank.getInt("sector", i);
            double layer     = adcBank.getInt("layer", i);
            double component = adcBank.getInt("component", i);
            double t0        = ahdcConstants.get_t0(sector, layer, component).t0;
            double leadingEdgeTime = adcBank.getFloat("leadingEdgeTime", i);
            //double timeMax         = adcBank.getFloat("time", i);
            double timeOverThreshold = adcBank.getFloat("timeOverThreshold", i);
            int    wfType  = adcBank.getInt("wfType", i);
            int    adcMax  = adcBank.getInt("ADC", i);
            int    adcOffset  = adcBank.getFloat("ped", i);
            double time       = leadingEdgeTime - t0;
            for (int type = 0; type <= 6; type++) {
                // wfType == type
                if (wfType == type) {
                    VecH1_wfType[type]->Fill(wfType); 
                    VecH1_time[type]->Fill(time); 
                    VecH1_leadingEdgeTime[type]->Fill(leadingEdgeTime); 
                    VecH1_timeOverThreshold[type]->Fill(timeOverThreshold); 
                    VecH1_adcMax[type]->Fill(adcMax); 
                    VecH1_adcOffset[type]->Fill(adcOffset); 
                    VecH1_occupancy[type]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component)); 
                }
                if (wfType <= type) {
                    VecH1_wfType[7+type]->Fill(wfType); 
                    VecH1_time[7+type]->Fill(time); 
                    VecH1_leadingEdgeTime[7+type]->Fill(leadingEdgeTime); 
                    VecH1_timeOverThreshold[7+type]->Fill(timeOverThreshold); 
                    VecH1_adcMax[7+type]->Fill(adcMax); 
                    VecH1_adcOffset[7+type]->Fill(adcOffset); 
                    VecH1_occupancy[7+type]->Fill(AhdcCCDB::wireUniqueId(sector, layer, component)); 
                }
            }
            if ((wfType < 0) || (wfType > 6)) printf("> WEIRD wfType : %d", wfType);
        }
    }
    printf("nevents    : %ld \n", nevents);
    
    // Renormalise occupancy
    for (int type = 0; type <= 6; type++) {
        VecH1_occupancy[type]->Scale(100.0/nevents); 
        VecH1_occupancy[7+type]->Scale(100.0/nevents); 
    }
    // output
    TFile *f = new TFile(output, "RECREATE");
   
    std::vector<TDirectory*> wfType_DIR;
    wfType_DIR.push_back(f->mkdir("wfType==0"));
    wfType_DIR.push_back(f->mkdir("wfType==1"));
    wfType_DIR.push_back(f->mkdir("wfType==2"));
    wfType_DIR.push_back(f->mkdir("wfType==3"));
    wfType_DIR.push_back(f->mkdir("wfType==4"));
    wfType_DIR.push_back(f->mkdir("wfType==5"));
    wfType_DIR.push_back(f->mkdir("wfType==6"));
    wfType_DIR.push_back(f->mkdir("wfType<=0"));
    wfType_DIR.push_back(f->mkdir("wfType<=1"));
    wfType_DIR.push_back(f->mkdir("wfType<=2"));
    wfType_DIR.push_back(f->mkdir("wfType<=3"));
    wfType_DIR.push_back(f->mkdir("wfType<=4"));
    wfType_DIR.push_back(f->mkdir("wfType<=5"));
    wfType_DIR.push_back(f->mkdir("wfType<=6"));

    // Loop over type
    for (int type = 0; type <= 6; type++) {
        wfType_DIR[type]->cd();
            VecH1_wfType[type]->Write(VecH1_wfType[type]->GetName());
            VecH1_time[type]->Write(VecH1_time[type]->GetName());
            VecH1_leadingEdgeTime[type]->Write(VecH1_leadingEdgeTime[type]->GetName());
            VecH1_timeOverThreshold[type]->Write(VecH1_timeOverThreshold[type]->GetName());
            VecH1_adcMax[type]->Write(VecH1_adcMax[type]->GetName());
            VecH1_adcOffset[type]->Write(VecH1_adcOffset[type]->GetName());
            VecH1_occupancy[type]->Write(VecH1_occupancy[type]->GetName());
        wfType_DIR[7+type]->cd(); 
            VecH1_wfType[7+type]->Write(VecH1_wfType[7+type]->GetName());
            VecH1_time[7+type]->Write(VecH1_time[7+type]->GetName());
            VecH1_leadingEdgeTime[7+type]->Write(VecH1_leadingEdgeTime[7+type]->GetName());
            VecH1_timeOverThreshold[7+type]->Write(VecH1_timeOverThreshold[7+type]->GetName());
            VecH1_adcMax[7+type]->Write(VecH1_adcMax[7+type]->GetName());
            VecH1_adcOffset[7+type]->Write(VecH1_adcOffset[7+type]->GetName());
            VecH1_occupancy[7+type]->Write(VecH1_occupancy[7+type]->GetName());
    }

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
