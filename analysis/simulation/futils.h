// futils : f stands for felix

#ifndef FUTILS_H
#define FUTILS_H

#include "dictionary.h"
#include <vector>

typedef double (*Integrable)(double);

/*
class Integrable {
    private :
        double (*ptr) (double) = nullptr;
        void get_ptr(double (*ptr_) (double)) const {ptr_=ptr; return;}
    public :
        Integrable(double (*ptr_) (double)){ptr = ptr_;}
        Integrable(const Integrable & func){func.get_ptr(ptr);}
        Integrable(){ptr=nullptr;}
        ~Integrable(){;}
        
        virtual double operator()(double x){
            return (*ptr)(x);
        }
};*/

namespace futils {

    // Statistical analysis
    double mean(const std::vector<double> & data);
    double variance(const std::vector<double> & data);
    double std_dev(const std::vector<double> & data);
    bool cart2polar(double x, double y, double z, double & rho, double & theta, double & phi);
    bool cart2polar(double x, double y, double & rho, double & theta);
    double integrate(Integrable f, double a, double b, int Npts);
}




#endif