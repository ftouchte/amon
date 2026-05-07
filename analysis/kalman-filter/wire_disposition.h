#ifndef WIRE_DISPOSITION_H
#define WIRE_DISPOSITION_H

#include <vector>

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

#endif