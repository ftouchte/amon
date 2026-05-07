/***********************************************
 * Look at 10.6 GeV data
 *
 * @author Felix Touchte Codjo
 * @date April 02, 2026
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

#include "lookat10gev.h"
#include "fOptions.h"
#include "Units.h"

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


int main(int argc, char const *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();
    
    /// --- Load files from options
    fOptions OPT({"-i", "-o"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");


    /// --- Histograms
    Histograms* Histos = new Histograms();


    /// --- Event counters
    long unsigned int nevents =0;

    /// --- Loop over files
    int nfiles = 0;
    for (std::string filename : filenames) {
        /// --- Open files
        nfiles++;
        printf("> Open filename %d/%d: %s\n", nfiles, (int) filenames.size(), filename.c_str());
        hipo::reader  reader(filename.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);
        
        /// --- Bank definition
        hipo::bank  adcBank(factory.getSchema("AHDC::adc"));
        hipo::bank  wfBank(factory.getSchema("AHDC::wf"));
        hipo::bank  trackBank(factory.getSchema("AHDC::kftrack"));
        hipo::bank  hitBank(factory.getSchema("AHDC::hits"));
        hipo::bank  recBank(factory.getSchema("REC::Particle"));
        hipo::event event;

        /// --- Loop over events
        while( reader.next()){
            
            /// --- display progress Bar
            nevents++;
            if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
                futils::progressBar(100.0*nevents/reader.getEntries());
            }

            /// --- load bank content for this event
            reader.read(event);
            event.getStructure(adcBank);
            event.getStructure(wfBank);
            event.getStructure(trackBank);
            event.getStructure(hitBank);
            event.getStructure(recBank);

            /// --- Loop over track
            for (int i = 0; i < trackBank.getRows(); i++) {
                int trackid = trackBank.getInt("trackid", i);
                int nhits = trackBank.getInt("n_hits", i);

                if (nhits < 6) continue;
                
                double vz = trackBank.getFloat("z", i) * Units::mm;
                double px = trackBank.getFloat("px", i) * Units::MeV;
                double py = trackBank.getFloat("py", i) * Units::MeV;
                double pz = trackBank.getFloat("pz", i) * Units::MeV;
                double dEdx = trackBank.getFloat("dEdx", i) * Units::MeV/Units::mm;
                // double dEdx = trackBank.getFloat("dEdx", i);

                double p = sqrt(px*px + py*py + pz*pz);
                double theta = acos(pz/p);
                double phi = atan2(py, px);

                Histos->H1_track_p->Fill(p);
                Histos->H1_track_theta->Fill(theta / Units::deg); // read per degree to convert in degree
                Histos->H1_track_phi->Fill(phi > 0 ? phi / Units::deg : 360 + phi / Units::deg);
                Histos->H1_track_nhits->Fill(nhits);
                Histos->H1_track_vz->Fill(vz);
                Histos->H2_track_corr_p_dEdx->Fill(p, dEdx / (Units::MeV/Units::mm));

                // Loop over hit for this track
                for (int j = 0; j < hitBank.getRows(); j++) {
                    if (hitBank.getInt("trackid", j) == trackid) {
                        double residual = hitBank.getDouble("residual", j) * Units::mm;
                        // double doca = hitBank.getDouble("doca", j) * Units::mm;
                        // double time = hitBank.getDouble("time", j) * Units::ns;
                        
                        Histos->H1_hit_residual->Fill(residual / Units::mm);
                    }
                }

            }

            /// --- Loop over electron
            for (int i = 0; i < recBank.getRows(); i++) {
                // Look for electrons
                if (recBank.getInt("pid", i) == 11) {
                    //int status = trackBank.getShort("status", i);
                    double vz = recBank.getFloat("vz", i) * Units::cm;
                    double px = recBank.getFloat("px", i) * Units::GeV;
                    double py = recBank.getFloat("py", i) * Units::GeV;
                    double pz = recBank.getFloat("pz", i) * Units::GeV;

                    double p = sqrt(px*px + py*py + pz*pz);
                    double theta = acos(pz/p);
                    double phi = atan2(py, px);

                    Histos->H1_electron_p->Fill(p);
                    Histos->H1_electron_theta->Fill(theta / Units::deg); // read per degree to convert in degree
                    Histos->H1_electron_phi->Fill(phi > 0 ? phi / Units::deg : 360 + phi / Units::deg);
                    Histos->H1_electron_vz->Fill(vz);
                    
                }
            }

            


        } // end loop over events for this file


    } // loop over files



    

    /// --- Output
    TFile *f = new TFile(output.c_str(), "RECREATE");
    Histos->SaveIn(f);


    f->Close();
    printf("> Output file created : %s\n", output.c_str());
    
    // end of the program
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("\033[1m * time elapsed : %lf seconds\033[0m\n", elapsed.count());
    return 0;
}


