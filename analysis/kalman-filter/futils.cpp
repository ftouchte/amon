// futils : f stands for felix

#define _USE_MATH_DEFINES

#include "./futils.h"
#include <iostream>
#include <cmath>


namespace futils {

    double mean(const std::vector<double> & data){
        long unsigned int N = data.size();
        double res=0;
        for (long unsigned int i=0; i<N;++i){
            res+= data.at(i);
        }
        return res/N;
    }

    double variance(const std::vector<double> & data){
        long unsigned int N = data.size();
        double res=0;
        double mean = futils::mean(data);
        for (long unsigned int i=0; i<N;++i){
            res+= pow(data.at(i) - mean,2);
        }
        return res/(N-1);
    }

    double std_dev(const std::vector<double> & data){
        return sqrt(futils::variance(data));
    }

    bool cart2polar(double x, double y, double z, double & rho, double & theta, double & phi){
        rho = sqrt(x*x+y*y+z*z);
        if (rho <= 1e-9) {theta = 0; phi = 0; return false;} // we consider that if x < 1e-9 then rho is 0
        if (sqrt(x*x+y*y) <= 1e-9) {phi = 0; return false;}
        theta = acos(z/rho);
        if (y >= 0){
            phi = acos(x/(rho*sin(theta)));
        }
        else {
            phi = 2*M_PI - acos(x/(rho*sin(theta)));
        }
        return true;
    } 

    bool cart2polarDEG(double x, double y, double z, double & rho, double & theta, double & phi){
        rho = sqrt(x*x+y*y+z*z);
        if (rho <= 1e-9) {theta = 0; phi = 0; return false;}
        if (sqrt(x*x+y*y) <= 1e-9) {phi = 0; return false;}
        theta = acos(z/rho);
        if (y >= 0){
            phi = acos(x/(rho*sin(theta)));
        }
        else {
            phi = 2*M_PI - acos(x/(rho*sin(theta)));
        }
        phi *= 180/M_PI;
        theta *= 180/M_PI;
        return true;
    } 

    bool cart2polar(double x, double y, double & rho, double & theta){
        rho = sqrt(x*x+y*y);
        if (rho <= 1e-9) {theta = 0; return false;}
        if (y >= 0) {
            theta = acos(x/rho);
        }
        else {
            theta = 2*M_PI - acos(x/rho);
        }
        return true;
    }

    bool cart2polarDEG(double x, double y, double & rho, double & theta){
        rho = sqrt(x*x+y*y);
        if (rho <= 1e-9) {theta = 0; return false;}
        if (y >= 0) {
            theta = acos(x/rho);
        }
        else {
            theta = 2*M_PI - acos(x/rho);
        }
        theta *= 180/M_PI;
        return true;
    }

    double integrate(Integrable f, double a, double b, int Npts){
        double S=0;
        double eps = (b-a)/Npts;
        for (int i=0;i<Npts;i++){
            S+= f(a+eps*i)+f(a+eps*(i+1)); // trapèze
        }
        return S*eps/2; // trapèze
    }
    

}
