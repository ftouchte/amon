/*******************************************************
 * fColorPalette.h
 *
 * @author Felix Touchte Codjo
 * @date February 17, 2025
 * ****************************************************/

#ifndef F_COLOR_PALETTE_H
#define F_COLOR_PALETTE_H

struct fColor {
	double r;
	double g;
	double b;
};

class fColorPalette {
	int n;
	std::vector<fColor> palette;
public :
	fColorPalette(int _n);
};


#endif
