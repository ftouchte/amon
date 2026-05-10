#ifndef WIRE_DISPOSITION_H
#define WIRE_DISPOSITION_H

#include <vector>
#include <map>

#include "TH1.h"
#include "TH2.h"

#include "Units.h"

struct Hit {
    int sector;
    int superlayer;
    int layer;
    int component;
    int wireNum;
    float adc;
    float time;
    float residual;
    float residual_LR;
};

struct Track {
    float px;
    float py;
    float pz;
    float p;
    float theta;
    float phi;
    int nhits;
    std::vector<Hit> hits;
    float electron_px;
    float electron_py;
    float electron_pz;
    float electron_p;
    float electron_theta;
    float electron_phi;
};


struct SingleThreadRun {

    /// --- Constants
    double W2_min = 3.5 * Units::GeV * Units::GeV;
    double W2_max = 3.8 * Units::GeV * Units::GeV;
    double vzInf_min = -20 * Units::cm;
    double vzInf_max = -10 * Units::cm;
    double vzSup_min = 5 * Units::cm;
    double vzSup_max = 15 * Units::cm;

    std::vector<int> all_bad_wires = {24, 25, 34, 35, 46, 47, 75, 76, 77, 88, 89, 90, 100, 102, 131, 132, 133, 140, 143, 144, 145, 146, 157, 165, 167, 196, 198, 214, 228, 369};
    std::vector<int> fixed_wires = {64, 65, 122, 124, 120};
    
    std::vector<int> bad_wires_vec1 = {24, 25, 75, 76, 77, 131, 132, 133};
    std::vector<int> bad_wires_vec2 = {34, 35, 88, 89, 90, 144, 145, 214};
    std::vector<int> bad_wires_vec3 = {47, 46, 102, 100, 226, 298, 157, 45, 158, 101};

    /// --- Stats
    int nfiles;
    long unsigned int nevents;
    long unsigned int nelectrons;

    /// --- Histograms
    TH1D* H1_W2;
    TH1D* H1_track_vz;

    std::vector<TH2D*> H2_corr_wire_phi;
    std::vector<TH2D*> H2_corr_wire_phi_vzInf;
    std::vector<TH2D*> H2_corr_wire_phi_vzSup;

    TH2D* H2_wire_occupancy;
    TH2D* H2_wire_occupancy_vzInf;
    TH2D* H2_wire_occupancy_vzSup;

    std::map<int,TH1D*> bad_wires_map1;
    std::map<int,TH1D*> bad_wires_map2;
    std::map<int,TH1D*> bad_wires_map3;

    std::map<int, TH1D*> leadingEdgeTimeBadWires_map;
    std::map<int,TH1D*>  residual_LRBadWires_map;

    std::vector<TH1D*> H1_all_residual_LR;
    std::vector<TH1D*> H1_all_leadingEdgeTime;

    ///< Constructor
    SingleThreadRun();

    ///< Destructor
    ~SingleThreadRun();

    void run(std::vector<std::string> files);

    void merge(const SingleThreadRun & t);

    void delete_map(std::map<int,TH1D*> & map) {
        for (const auto& [wireNum, h] : map) {
            delete h;
        }
        map.clear();
    }

    template<typename T>
    void delete_histo_vector(std::vector<T*> & vec) {
        for (auto h : vec) {
            delete h;
        }
        //vec.clear();
    }

private:
    void add_map_histo(std::map<int,TH1D*> & map_dest, const std::map<int,TH1D*> & map_src);

    template<typename T>
    void add_histo_vector(std::vector<T*> & vec_dest, const std::vector<T*> & vec_src) {
        for (size_t i = 0; i < vec_dest.size(); i++) {
            vec_dest[i]->Add(vec_src[i]);
        }
    }
    
};

#endif