/***********************************************
 * test_Histograms.cpp
 *
 * @author Felix Touchte Codjo
 * @date February 13, 2025
 * ********************************************/

#include <cstdio>
#include <random>
#include "fH1D.h"


int main(int argc, const char * argv[]) {
	printf("=====> Test Histogram 1D\n");		
	// gauss
	fH1D hist1d_gauss("gauss", 50, -3, 3);
	fH1D hist1d_cauchy("cauchy", 50, -3,3);
	fH1D hist1d_lognormal("log normal", 50, -1, 5);
	std::random_device rd{};
	std::mt19937 gen{rd()};
	std::normal_distribution dg{0.0, 1.0};
	std::cauchy_distribution<float> dc{0.0, 1.0};
	std::lognormal_distribution<> dlg(0.0, 1.5);
	for (int i = 0; i < 1'000'000; i++) {
		hist1d_gauss.fill(dg(gen));
		hist1d_cauchy.fill(dc(gen));
		hist1d_lognormal.fill(dlg(gen));
	}
	hist1d_gauss.print();
	hist1d_cauchy.print();
	hist1d_lognormal.print();
		
	return 0;
}
