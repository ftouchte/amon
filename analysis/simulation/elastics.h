#ifndef ELASTIC_H
#define ELASTIC_H

struct State {
    double px, py, pz;
    double vx, vy, vz;
    double p, theta, phi;
    double pT;
    int pid = -9999; // only of rec::particle
    State(double _px, double _py, double _pz, double _vx, double _vy, double _vz);
    // only relavant for track
    int trackid = -9999;
    double dEdx = -9999;
    double adc = -9999;
};

struct Physics {
    double Q2;
    double W2;
    double nu;
    double xB;
    Physics(State Probe, double _Mt, double _Ee, double _me = 0.511e-3, double _Mp = 938.272e-3);
};

struct ElasticsOutput {
    State probe;
    State track;
};

#endif
