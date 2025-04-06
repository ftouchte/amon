/**************************************
 * f3Dutils.h
 *
 * @author Felix Touchte Codjo
 * @date April 6, 2025
 * ***********************************/


#ifndef F_3D_UTILS_H
#define F_3D_UTILS_H

#include <gtkmm.h>

struct V3D {
	double x;
	double y;
	double z;
};

struct PseudoV3D {
	double X;
	double Y;
};

class Frame {
	double alpha;
	double beta;
	double gamma;
	PseudoV3D vx;
	PseudoV3D vy;
	PseudoV3D vz;
public :
	Frame(double _alpha, double _beta, double _gamma);

	PseudoV3D getXAxis();
	PseudoV3D getYAxis();
	PseudoV3D getZAxis();

	PseudoV3D get_PseudoV3D(double x, double y, double z);
	void update();
	void set_alpha(double _alpha);
	void set_beta (double _beta);
	void set_gamma(double _gamma);

	// test
	static void test1(double alpha, double beta, double gamma, const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);	
};


#endif
