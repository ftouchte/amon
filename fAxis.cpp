/************************************************* 
 * Class for axis (labels and sticks)
 * 
 * @author Felix Touchte Codjo 
 * @date January 26,2025 
 * **********************************************/

#include "fAxis.h"
#include <cmath>
#include <cstdio>


// Utility
// int sprintf( char* buffer, const char* format, ... );

/** Constructor */
fAxis::fAxis(double _start, double _end, int _n1, int _n2, int _n3) : start(_start), end(_end), n1(_n1), n2(_n2), n3(_n3) {
	width = end - start;
	if (width < 0) { n1 = 0;}
	/********************
	 * first divisions
	 * ******************/
	if (n1 > 0) {
		eps1 = width/n1;
		// what is closer ?
		// - 1*10^alpha 
		// - 2*10^beta
		// - 5*10^gamma
		int alpha = ceil(log10(eps1/1.0));
		int beta  = ceil(log10(eps1/2.0));
		int gamma = ceil(log10(eps1/5.0));
		std::vector<double> precisions = {1*pow(10.0,alpha), 2*pow(10.0,beta), 5*pow(10.0,gamma)};
		std::vector<int> tab_125 = {1, 2, 5};
		std::vector<int> powers = {alpha, beta, gamma};
		double closer = precisions[0] - eps1;
		int icloser = 0;
		for (int i = 0; i < (int) precisions.size(); i++) {
			if (closer > (precisions[i] - eps1)) {
				closer  = precisions[i] - eps1;
				icloser = i;
			}
		}
		//printf("(alpha, beta, gamma) : (%d, %d, %d)\n", alpha, beta, gamma);
		//printf("precision : %d*10^%d\n", tab_125[icloser], powers[icloser]);
		// we redefine eps1 
		// n_inf*eps1 < start < end < n_sup*eps1
		eps1 = tab_125[icloser]*pow(10.0, powers[icloser]);
		int ndecimals = powers[icloser] < 0 ? -powers[icloser] : 0; 
		int n_inf = floor(start/eps1);
		int n_sup = ceil(end/eps1);
		// list of first divisions
		for (int k = n_inf; k <= n_sup; k++) {
			double value = k*eps1;
			Div1.push_back(value);
			char buffer[50];
			sprintf(buffer, "%.*lf", ndecimals, value);
			Labels1.push_back(buffer);
		}
		/**********************
		 * second divisions
		 * ********************/
		if (n2 > 0) {
			double eps2 = eps1/n2;
			for (int k = 1; k <= (int) Div1.size() - 1; k++) {
				for (int i = 1; i <= n2-1; i++) {
					double  value = Div1[k-1] + i*eps2;
					Div2.push_back(value);
					char buffer[50];
					sprintf(buffer, "%.*lf", (value == floor(value) ? 0 : ndecimals + 1), value);
					Labels2.push_back(buffer);
				}
			}
			/**********************
			 * second divisions
			 * ********************/
			if (n3 > 0) {
				eps3 = eps2/n3;
				for (int k = 1; k <= (int) Div2.size() - 1; k++) {
					for (int i = 1; i < n3; i++) {
						double value = Div2[k-1] + i*eps3;
						Div2.push_back(value);
						char buffer[50];
						sprintf(buffer, "%.*lf", (value == floor(value) ? 0 : ndecimals + 2), value);
						Labels3.push_back(buffer);
					}
				}
			}
		}

	}
}


fAxis::~fAxis(){} ///< Destructor

int  fAxis::get_n1() {return n1;} ///< get the number of first divisions
int  fAxis::get_n2() {return n2;} ///< get the number of second divisions
int  fAxis::get_n3() {return n3;} ///< get the number of third divisions
int  fAxis::get_eps1() {return eps1;} ///< get the space between the first divisions
int  fAxis::get_eps2() {return eps2;} ///< get the space between the second divisions
int  fAxis::get_eps3() {return eps3;} ///< get the space between the third divisions

void fAxis::print() {
	printf(">>> start : %lf, end : %lf, eps1 : %lf\n",start,end,eps1);
	printf("    1st divisions [ ");
	for (std::string s : Labels1) {
		printf("%s ", s.c_str());
	}
	printf("]\n");

	printf("    2nd divisions [ ");
	for (std::string s : Labels2) {
		printf("%s ", s.c_str());
	}
	printf("]\n");

	printf("    3rd divisions [ ");
	for (std::string s : Labels3) {
		printf("%s ", s.c_str());
	}
	printf("]\n");
}


std::vector<std::string> fAxis::get_labels1() {
	return Labels1;
}


std::vector<std::string> fAxis::get_labels2() {
	return Labels2;
}

std::vector<std::string> fAxis::get_labels3() {
	return Labels3;
}
