/**********************************************
 * Lund file generator
 *
 * cf. https://gemc.jlab.org/gemc/html/documentation/generator/lund.html
 *
 * @author : ftouchte
 * @date : November 25, 2025
 * *******************************************/

#include <fstream>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <random>

#include "TH1.h"
#include "TFile.h"
#include "TString.h"

struct LundHeader {
    //***  number of particle to generate in this event
    int nbOfParticles;
    // mass number of the target (UD)
    int A;
    // atomic number (UD)
    int Z;
    // target polarisation (UD)
    double polarization = 0; 
    //*** first particle Z component of spin
    double spin = 0; 
    // 11 for electron and 22 for photon (UD)
    int beam_type = 11;
    // beam energy (GeV) (UD)
    double beam_energy; 
    // interacted nucleon ID (2212 or 2112) (UD)
    int nucleon_id = 2212;
    // process ID (UD)
    int process_id = 1;
    // event weight (UD)
    double event_weight = 1;

};

struct LundParticle {
    int index;
    double lifetime = -1; // UD
    int type = 1;
    int ID;
    int parent_index = 0;
    int daughter_index = 0; // UD
    // momentum [GeV]
    double px;
    double py;
    double pz;
    double p;
    double theta; // radians
    double phi; // radians
    double energy = 0; // UD
    double mass = 0; // UD
    // vertex [cm]
    double vx;
    double vy;
    double vz;
};

void polar2cart(double p, double theta, double phi, double & px, double & py, double & pz) {
    pz = p*cos(theta);
    px = p*sin(theta)*cos(phi);
    py = p*sin(theta)*sin(phi);
}

double uniform(double a, double b, int nevents) {
    double x = (1.0*(rand() % nevents))/nevents;
    return (b-a)*x + a;
}

void generate_kinematics(int nevents, LundParticle & electron, LundParticle & track) {
    // beam energy
    double Ee = 2.23951; // GeV
    double nu = 0.06966; // GeV, nu = E - E', this the mean of the distribution
    // electron
    double phi = uniform(0,360, nevents)*M_PI/180; // deg to rad
    double theta = uniform(7,8, nevents)*M_PI/180;
    double p = Ee - nu; // GeV
    
    polar2cart(p, theta, phi, electron.px, electron.py, electron.pz);
    electron.p = p;
    electron.theta = theta;
    electron.phi = phi;
    // track
    double track_p = 2*Ee*sin(theta/2);
    double track_pz = Ee*(1-cos(theta));
    double track_theta = acos(track_pz/track_p);
    double track_phi = (phi > M_PI) ? (phi - M_PI) : (phi + M_PI);
    polar2cart(track_p, track_theta, track_phi, track.px, track.py, track.pz);
    track.p = track_p;
    track.theta = track_theta;
    track.phi = track_phi;
    // vertex
    double vertex = uniform(-15,15, nevents); // cm
    electron.vx = 0;
    electron.vy = 0;
    electron.vz = vertex;
    track.vx = 0;
    track.vy = 0;
    track.vz = vertex;
}

#include <filesystem>
int main(int argc, const char * argv[]) {
    
    //const char * output_name = "lund_tobesimulated.dat";
    //std::ofstream ofs (output_name, std::ofstream::out);
    if (argc < 1) {
        printf("Please, provide an empty directory...\n");
        return 0;
    }
    const char * output_dir = argv[1];
    if (!(std::filesystem::exists(output_dir) && std::filesystem::is_directory(output_dir) && std::filesystem::is_empty(output_dir))) {
        printf("Error: (%s) path does not exist or path is not a directory or directory is not empty.\n", output_dir);
        return 0;
    }
    // total number of events = nevents*nfiles
    int nevents = 2000;
    int nfiles = 40;
    std::vector<std::ofstream> ofs;
    for (int i = 0; i < nfiles; i++) {
        ofs.push_back(std::ofstream(TString::Format("%s/lund_file_%0*d.dat", output_dir,(int) floor(log10(nfiles))+1, i).Data(), std::ofstream::out));
    }
    
    // Header
    LundHeader header;
    header.nbOfParticles = 2;
    header.A = 2;
    header.Z = 1;
    header.polarization = 0; 
    header.spin = 0; 
    header.beam_type = 11;
    header.beam_energy = 2.24; 
    header.nucleon_id = 2212;
    header.process_id = 1;
    header.event_weight = 1;
    // electron
    LundParticle electron;
    electron.index = 1;
    electron.lifetime = -1; // UD
    electron.type = 1;
    electron.ID = 11;
    electron.parent_index = 0;
    electron.daughter_index = 0; // UD
    electron.px = 0.24;
    electron.py = 0;
    electron.pz = 2.17;
    electron.mass = 0.511e-3; // UD
    electron.energy = sqrt(pow(electron.mass,2) + pow(electron.px,2) + pow(electron.py,2) + pow(electron.pz,2)); // UD
    electron.vx = 0;
    electron.vy = 0;
    electron.vz = 0;
    // track
    LundParticle track;
    track.index = 2;
    track.lifetime = -1; // UD
    track.type = 1;
    track.ID = 1000010020; //   this is a deuteron
    track.parent_index = 0;
    track.daughter_index = 0; // UD
    track.px = -0.24;
    track.py = 0;
    track.pz = 0.03;
    track.mass = 1.875;
    track.energy = sqrt(pow(track.mass,2) + pow(track.px,2) + pow(track.py,2) + pow(track.pz,2)); // UD
    track.vx = 0;
    track.vy = 0;
    track.vz = 0;
    
    // Histograms
    TH1D* H1_electron_p = new TH1D("electron_p", "p; p (GeV); count", 100, 2.17, 2.20);
    TH1D* H1_electron_theta = new TH1D("electron_theta", "theta; theta (deg); count", 100, 0, 90);
    TH1D* H1_electron_phi = new TH1D("electron_phi", "phi; phi (deg); count", 100, 0, 361);
    TH1D* H1_electron_vz = new TH1D("electron_vz", "vz; vz (cm); count", 100, -15, 15);
    TH1D* H1_track_p = new TH1D("track_p", "p; p (GeV); count", 100, 0.15, 0.5);
    TH1D* H1_track_theta = new TH1D("track_theta", "theta; theta (deg); count", 100, 0, 90);
    TH1D* H1_track_phi = new TH1D("track_phi", "phi; phi (deg); count", 100, 0, 361);
    TH1D* H1_track_vz = new TH1D("track_vz", "vz; vz (cm); count", 100, -15, 15);
    // Loop over portions of file
    for (int num = 0; num < nfiles; num++) {
        // Loop over number of events
        for (int i = 0; i <= nevents; i++) {
            // Update kinematics
            generate_kinematics(nfiles*nevents, electron, track);
            H1_electron_p->Fill(electron.p);
            H1_electron_theta->Fill(electron.theta*180/M_PI);
            H1_electron_phi->Fill(electron.phi*180/M_PI);
            H1_electron_vz->Fill(electron.vz);
            H1_track_p->Fill(track.p);
            H1_track_theta->Fill(track.theta*180/M_PI);
            H1_track_phi->Fill(track.phi*180/M_PI);
            H1_track_vz->Fill(track.vz);
            // Ici l'ordre est très important !!!
            // Print Header
            ofs[num] << header.nbOfParticles << " ";
            ofs[num] << header.A << " ";
            ofs[num] << header.Z << " ";
            ofs[num] << header.polarization << " "; 
            ofs[num] << header.spin << " "; 
            ofs[num] << header.beam_type << " ";
            ofs[num] << header.beam_energy << " "; 
            ofs[num] << header.nucleon_id << " ";
            ofs[num] << header.process_id << " ";
            ofs[num] << header.event_weight << std::endl;
            // Print electron
            ofs[num] << electron.index << " ";
            ofs[num] << electron.lifetime << " "; // UD
            ofs[num] << electron.type << " ";
            ofs[num] << std::setw(10) << electron.ID << " ";
            ofs[num] << electron.parent_index << " ";
            ofs[num] << electron.daughter_index << " "; // UD
            ofs[num] << electron.px << " ";
            ofs[num] << electron.py << " ";
            ofs[num] << electron.pz << " ";
            ofs[num] << electron.energy << " "; // UD
            ofs[num] << electron.mass << " "; // UD
            ofs[num] << electron.vx << " ";
            ofs[num] << electron.vy << " ";
            ofs[num] << electron.vz << std::endl;
            // Print track
            ofs[num] << track.index << " ";
            ofs[num] << track.lifetime << " "; // UD
            ofs[num] << track.type << " ";
            ofs[num] << std::setw(10) << track.ID << " ";
            ofs[num] << track.parent_index << " ";
            ofs[num] << track.daughter_index << " "; // UD
            ofs[num] << track.px << " ";
            ofs[num] << track.py << " ";
            ofs[num] << track.pz << " ";
            ofs[num] << track.energy << " "; // UD
            ofs[num] << track.mass << " "; // UD
            ofs[num] << track.vx << " ";
            ofs[num] << track.vy << " ";
            ofs[num] << track.vz << std::endl;
        }
    } 
    // save in root file
    TFile *f = new TFile("./output/lund_file.root", "RECREATE");
    H1_electron_p->Write("electron_p");
    H1_electron_theta->Write("electron_theta");
    H1_electron_phi->Write("electron_phi");
    H1_electron_vz->Write("electron_vz");
    H1_track_p->Write("track_p");
    H1_track_theta->Write("track_theta");
    H1_track_phi->Write("track_phi");
    H1_track_vz->Write("track_vz");
    f->Close();
    printf("File created : %s\n", "./output/lund_file.root");
    
    
    // close file
    for (int i = 0; i < nfiles; i++) {
        ofs[i].close();
        printf("File created : %s\n", TString::Format("./output-lund/lund_file_%0*d.dat", (int) floor(log10(nfiles))+1, i).Data());
    }

    return 0;

}
