/**********************************************
 * Class to laod a ahdc ccdb constant
 *
 * @author ftouchte
 * @date August 17, 2025
 * ********************************************/

#include "AhdcCCDB.h"

#include <string>
#include <cstdio>
#include <vector>

#include <CCDB/Calibration.h>
#include <CCDB/CalibrationGenerator.h>
#include <CCDB/SQLiteCalibration.h>

AhdcCCDB::AhdcCCDB(std::string _connection, int _runNo, std::string _variation, std::string _timestamp) 
    : connection(_connection), runNo(_runNo), variation(_variation), timestamp(_timestamp) {
    
    // Create calibraion object
    // To specify run number, variation and time, the full signature of CreateCalibration is:
    // CreateCalibration(connection_str, int run=0, const string& variation="default", const time_t time=0);
    //
    // (!) remark: CreateCalibration throws logic_error with explanation if something goes wrong
    unique_ptr<ccdb::Calibration> calib(ccdb::CalibrationGenerator::CreateCalibration(connection));

    printf("********************** AHDC CCDB *****************************\n");
    printf("* > connection : %s\n", connection.c_str());
    printf("* > run number : %d\n", runNo);
    printf("* > variation  : %s\n", variation.c_str());
    printf("* > timestamp  : %s\n", timestamp.c_str());
    // The calibration object is now connected and ready to work. Lets get data:
    std::string dbname, database;
    std::vector<std::vector<double>> data;

    /////////////////////
    // raw_hit_cuts
    // //////////////////
    dbname = "/calibration/alert/ahdc/raw_hit_cuts";
    database = dbname + ":" + to_string(runNo) + ":" + variation + ":" + timestamp;
    calib->GetCalib(data, database.c_str());
    //iterate data and print
    for(int row = 0; row < (int) data.size(); row++) {
        int sector    = data[row][0];
        int layer     = data[row][1];
        int component = data[row][2];
        Cuts[wireUniqueId(sector, layer, component)] = {data[row][3], data[row][4], data[row][5], data[row][6], data[row][7], data[row][8], data[row][9], data[row][10]};
    }
    data.clear();
    printf("* > read %s\n", dbname.c_str());
    /////////////////////
    // time_offsets : t0
    // //////////////////
    dbname = "/calibration/alert/ahdc/time_offsets";
    database = dbname + ":" + to_string(runNo) + ":" + variation + ":" + timestamp;
    calib->GetCalib(data, database.c_str());
    //iterate data and print
    for(int row = 0; row < (int) data.size(); row++) {
        int sector    = data[row][0];
        int layer     = data[row][1];
        int component = data[row][2];
        T0Correction[wireUniqueId(sector, layer, component)] = {data[row][3], data[row][4], data[row][7]};
    }
    data.clear();
    printf("* > read %s\n", dbname.c_str());
    /////////////////////
    // time_to_distance
    // //////////////////
    dbname = "/calibration/alert/ahdc/time_to_distance";
    database = dbname + ":" + to_string(runNo) + ":" + variation + ":" + timestamp;
    calib->GetCalib(data, database.c_str());
    //iterate data and print
    for(int row = 0; row < (int) data.size(); row++) {
        T2d = {data[row][3], data[row][4], data[row][5], data[row][6], data[row][7], data[row][8],
               data[row][9], data[row][10], data[row][11], data[row][12], data[row][13], data[row][14],
               data[row][15],};
    }
    data.clear();
    printf("* > read %s\n", dbname.c_str());
    printf("********************** AHDC CCDB (end) ***********************\n");
}

int AhdcCCDB::wireUniqueId(int sector, int layer, int component) {
    if      (layer == 11) {
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

ahdcT0 AhdcCCDB::get_t0(int sector, int layer, int component) {
    return T0Correction[wireUniqueId(sector, layer, component)];
}

ahdcT0 AhdcCCDB::get_t0(int wire) {
    return T0Correction[wire];
}

ahdcRawCuts AhdcCCDB::get_rawCuts(int sector, int layer, int component) {
    return Cuts[wireUniqueId(sector, layer, component)];
}

ahdcRawCuts AhdcCCDB::get_rawCuts(int wire) {
    return Cuts[wire];
}

ahdcT2d AhdcCCDB::get_t2d() {
    return T2d;
}
