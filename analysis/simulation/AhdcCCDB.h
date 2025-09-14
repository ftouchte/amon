/**********************************************
 * Class to laod a ahdc ccdb constant
 *
 * @author ftouchte
 * @date August 17, 2025
 * ********************************************/


#ifndef AHDC_CCDB_H
#define AHDC_CCDB_H

#include <string>

struct ahdcRawCuts {
    double t_min;
    double t_max;
    double tot_min;
    double tot_max;
    double adc_min;
    double adc_max;
    double ped_min;
    double ped_max;
};

struct ahdcT2d {
    double p0;
    double p1;
    double p2;
    double p3;
    double p4;
    double p5;
    double dp0;
    double dp1;
    double dp2;
    double dp3;
    double dp4;
    double dp5;
    double chi2ndf;
};

struct ahdcT0 {
    double t0;
    double dt0;
    double chi2ndf;
};

class AhdcCCDB {
    std::string connection;
    int runNo;
    std::string variation;
    std::string timestamp; //e.g "2025-05-15_10-22-27"
   
    // t0 table
    ahdcT0 T0Correction[576];
    // raw cut tables
    ahdcRawCuts Cuts[576];
    // time to distance
    ahdcT2d T2d;
public:
    AhdcCCDB(std::string _connection = "mysql://clas12reader@clasdb.jlab.org/clas12", int _runNo = 22000, std::string _variation = "default", std::string timestamp = "no"); 
    static int wireUniqueId(int sector, int layer, int component);
    ahdcT0 get_t0(int sector, int layer, int component);
    ahdcT0 get_t0(int wire);
    ahdcRawCuts get_rawCuts(int sector, int layer, int component);
    ahdcRawCuts get_rawCuts(int wire);
    ahdcT2d get_t2d();
};

#endif
