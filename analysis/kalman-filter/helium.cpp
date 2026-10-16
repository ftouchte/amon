/***********************************************
 * Study elastics events off He4
 *
 * @author Felix Touchte Codjo
 * @date October 15, 2026
 * 
 * ********************************************/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vector>
#include <string>
#include <chrono>

#include "helium.h"
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
#include "TTree.h"
#include "TMath.h"

// nuclear masses
const double electron_mass = 0.511*Units::MeV; ///< electron energy mass
const double proton_mass = 938.272*Units::MeV; ///< proton energy mass
const double neutron_mass = 939.565*Units::MeV; ///< neutron energy mass
const double nucleon_mass = 0.5*(proton_mass + neutron_mass); ///< nucleon mass
const double helium_mass = 3.73*Units::GeV; ///< heluim4 energy mass
const double deuteron_mass = 1.875*Units::GeV; ///< deuteron energy mass

// beam energy
const double beam_energy =  2239.49*Units::MeV;


double mod2Pi(double x) {
    double r = atan2(sin(x), cos(x)); // is an angle between -pi and pi
    if (r >= 0) {
        return r;
    } else {
        return r + 2*M_PI;
    }
}

int main(int argc, char const *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();
    
    /// --- Load files from options
    fOptions OPT({"-i", "-o"});
    OPT.LoadOptions(argc, argv);
    OPT.Show();

    std::vector<std::string> filenames = OPT.GetValues("-i");
    std::string output = OPT.GetValue("-o");

    /// --- Histograms
    Histograms* histos = new Histograms();


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
        
        // Event counters
        long unsigned int nevents =0;

        /// --- Loop over events
        while( reader.next()){
            
            /// --- display progress Bar
            nevents++;
            if ((nevents % 1000 == 0) || ((int) nevents == reader.getEntries())) {
                progressBar(100.0*nevents/reader.getEntries());
            }

            /// --- load bank content for this event
            reader.read(event);
            event.getStructure(adcBank);
            event.getStructure(wfBank);
            event.getStructure(trackBank);
            event.getStructure(hitBank);
            event.getStructure(recBank);

            // kinematics kine;
            particle incoming_electron;
            {
                double E = beam_energy;
                double px = 0;
                double py = 0;
                double pz = sqrt(E*E - electron_mass*electron_mass);
                incoming_electron.momentum.SetCoordinates(px, py, pz, E);
            }
            particle initial_helium;
            {
                double E = helium_mass;
                double px = 0;
                double py = 0;
                double pz = 0;
                initial_helium.momentum.SetCoordinates(px, py, pz, E);
            }

            particle outgoing_electron;
            particle recoil_helium;
            particle virtual_photon;
            particle sys_eX;


            // Find one electron
            for (int row = 0; row < recBank.getRows(); row++) {
                if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0) {
                    // fill histograms
                    double px = recBank.getDouble("px", row)*Units::GeV;
                    double py = recBank.getDouble("py", row)*Units::GeV;
                    double pz = recBank.getDouble("pz", row)*Units::GeV;
                    double vz = recBank.getDouble("vz", row)*Units::cm;

                    // compute spherical coordiantes
                    double p = sqrt(px*px+py*py+pz*pz);
                    // double phi = mod2Pi(atan2(py, px));
                    // double theta = acos(pz/p);
                    double E = sqrt(p*p + electron_mass*electron_mass);

                    // // about the physics
                    double scattered_beam_energy = sqrt(pow(p,2) + pow(electron_mass,2));
                    double nu = beam_energy - scattered_beam_energy;
                    double y = nu/beam_energy;
                    // double Q2 = 4*beam_energy*scattered_beam_energy*pow(sin(theta/2),2);
                    // double W2 = pow(nucleon_mass,2) + 2*nucleon_mass*nu - Q2;
                    // double xB = Q2/(2*proton_mass*nu);

                    // // fill data and histo
                    // kine.Q2 = Q2;
                    // kine.xB = xB;
                    // kine.W2 = W2;
                    // kine.y = y;
                    // kine.nu = nu;

                    outgoing_electron.momentum.SetCoordinates(px, py, pz, E);
                    virtual_photon.momentum = incoming_electron.momentum - outgoing_electron.momentum;
                    sys_eX.momentum = virtual_photon.momentum + initial_helium.momentum;

                    double Q2 = -virtual_photon.momentum.M2();
                    double W2 = sys_eX.momentum.M2();
                    double xB = Q2/(2*initial_helium.momentum.Dot(virtual_photon.momentum));


                    // Fill histograms
                    histos->h1_Q2->Fill(Q2);
                    histos->h1_W2->Fill(W2);
                    histos->h1_xB->Fill(xB);
                    histos->h1_nu->Fill(nu);
                    histos->h1_y->Fill(y);
                    histos->h1_electron_vz_nocut->Fill(vz);

                    // Select eleastic based on the W2 cut
                    if (W2 > 13.8*Units::GeV2 && W2 < 14.04*Units::GeV2 && vz > -29 && vz < 13) {

                        outgoing_electron.row = row;

                        histos->h1_Q2_W2cut->Fill(Q2);
                        histos->h1_W2_W2cut->Fill(W2);
                        histos->h1_xB_W2cut->Fill(xB);
                        histos->h1_nu_W2cut->Fill(nu);
                        histos->h1_y_W2cut->Fill( y);

                        ThreeVector vec = outgoing_electron.momentum.Vect();
                        histos->h1_electron_p->Fill(vec.R());
                        histos->h1_electron_theta->Fill(vec.Theta()/Units::deg);
                        histos->h1_electron_phi->Fill(mod2Pi(vec.Phi())/Units::deg);
                        histos->h1_electron_vz->Fill(vz);
                    }


                    break;

                } // end pid == 1
            }

            // Find one track
            if (outgoing_electron.row < 0) continue;
            for (int row = 0; row < trackBank.getRows(); row++) {
                        
                int nhits = trackBank.getInt("n_hits", row);
                double chi2 = trackBank.getFloat("chi2", row);

                histos->h1_ahdc_track_nhits_nocut->Fill(nhits);
                histos->h1_ahdc_track_chi2_nocut->Fill(chi2);

                if (nhits >= 7 && chi2 < 8) {
                
                    double vz = trackBank.getFloat("z",  row) * Units::mm;
                    double px = trackBank.getFloat("px", row) * Units::MeV;
                    double py = trackBank.getFloat("py", row) * Units::MeV;
                    double pz = trackBank.getFloat("pz", row) * Units::MeV;
                    double dEdx = trackBank.getFloat("dEdx", row);

                    double p = sqrt(px*px + py*py + pz*pz);
                    double phi = mod2Pi(atan2(py, px));
                    double theta = acos(pz/p);
                    double E = sqrt(p*p + helium_mass*helium_mass);

                    histos->h1_ahdc_track_nhits->Fill(nhits);
                    histos->h1_ahdc_track_chi2->Fill(chi2);

                    histos->h1_ahdc_track_p->Fill(p);
                    histos->h1_ahdc_track_theta->Fill(theta/Units::deg);
                    histos->h1_ahdc_track_phi->Fill(mod2Pi(phi)/Units::deg);
                    histos->h1_ahdc_track_vz->Fill(vz);

                    recoil_helium.momentum.SetCoordinates(px, py, pz, E);

                    // Look at delta phi
                    ThreeVector v1 = outgoing_electron.momentum.Vect();
                    ThreeVector v2 = recoil_helium.momentum.Vect();
                    double phi1 = v1.Phi();
                    double phi2 = v2.Phi();
                    double delta_phi = mod2Pi(phi1- phi2) - M_PI;

                    histos->h1_ahdc_track_delta_phi->Fill(delta_phi/Units::deg);

                    if (delta_phi > -5*Units::deg && delta_phi < 5*Units::deg) {
                        histos->h1_elastics_ahdc_track_delta_phi->Fill(delta_phi/Units::deg);
                        histos->h1_elastics_ahdc_track_nhits->Fill(nhits);
                        histos->h1_elastics_ahdc_track_chi2->Fill(chi2);
                        histos->h1_elastics_ahdc_track_p->Fill(p);
                        histos->h1_elastics_ahdc_track_theta->Fill(theta/Units::deg);
                        histos->h1_elastics_ahdc_track_phi->Fill(mod2Pi(phi)/Units::deg);
                        histos->h1_elastics_ahdc_track_vz->Fill(vz);
                        histos->h2_elastics_ahdc_track_p_vs_dEdx->Fill(p, dEdx);
                        histos->h1_elastics_electron_p->Fill(v1.R());
                        histos->h1_elastics_electron_theta->Fill(v1.Theta()/Units::deg);
                        histos->h1_elastics_electron_phi->Fill(mod2Pi(v1.Phi())/Units::deg);
                        histos->h1_elastics_electron_vz->Fill(recBank.getDouble("vz", outgoing_electron.row));
                    }

                } // end cut track quality

            } // end loop over AHDC track






        } // end loop over events

    } // end loop over files













    // Output
    TFile *f = new TFile(output.c_str(), "RECREATE");

    histos->Write(f);

    f->Close();
    printf("> Output file created : %s\n", output.c_str());
    
    // End of the program
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start);
    printf("\033[1m * time elapsed : %lf seconds\033[0m\n", elapsed.count());
    return 0;

}
