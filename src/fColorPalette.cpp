/*******************************************************
 * fColorPalette.cpp
 *
 * @author Felix Touchte Codjo
 * @date March 16, 2025
 * ****************************************************/

#include "fColorPalette.h"

fColorPalette::fColorPalette(int _ndiv, int _type) : ndiv(_ndiv), type(_type) {
	if (type == 1) {
		int nrange = PALETTE1.size() - 1;
		for (int i = 1; i <= nrange; i++) {
			for (int j = 0; j < ndiv; j++) {
				double t = (1.0*j)/ndiv;	
				double r = (1-t)*PALETTE1[i-1].r + t*PALETTE1[i].r;
				double g = (1-t)*PALETTE1[i-1].g + t*PALETTE1[i].g;
				double b = (1-t)*PALETTE1[i-1].b + t*PALETTE1[i].b;
				color_range.push_back({r, g, b});
			}
		}
		color_range.push_back({PALETTE1[nrange].r, PALETTE1[nrange].g, PALETTE1[nrange].b});
	}
	else if (type == 2) {
		int nrange = PALETTE2.size() - 1;
		for (int i = 1; i <= nrange; i++) {
			for (int j = 0; j < ndiv; j++) {
				double t = (1.0*j)/ndiv;	
				double r = (1-t)*PALETTE2[i-1].r + t*PALETTE2[i].r;
				double g = (1-t)*PALETTE2[i-1].g + t*PALETTE2[i].g;
				double b = (1-t)*PALETTE2[i-1].b + t*PALETTE2[i].b;
				color_range.push_back({r, g, b});
			}
		}
		color_range.push_back({PALETTE2[nrange].r, PALETTE2[nrange].g, PALETTE2[nrange].b});
	}
	else {
		// do nothing
	}
}

fColor fColorPalette::get_color(int i) { return color_range.at(i);}
int fColorPalette::get_ncolors() { return color_range.size();}
std::vector<fColor> fColorPalette::get_color_range() { return color_range;}


