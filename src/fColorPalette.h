/*******************************************************
 * fColorPalette.h
 *
 * @author Felix Touchte Codjo
 * @date March 16, 2025
 * ****************************************************/

#ifndef F_COLOR_PALETTE_H
#define F_COLOR_PALETTE_H

#include <vector>

struct fColor {
	double r;
	double g;
	double b;
};

inline const std::vector<fColor> PALETTE1 = {
	// {r, g, b}
	{0.0, 0.0, 0.6}, // dark blue
	{0.0, 0.0, 1.0}, // blue
	{0.0, 1.0, 1.0}, // blue-sky
	{0.0, 1.0, 0.0}, // green
	{1.0, 1.0, 0.0}, // yellow
	{1.0, 0.0, 0.0}, // red
	{0.8, 0.0, 0.0} // dark red
};


inline const std::vector<fColor> PALETTE2 = {
	// {r, g, b}
	{0.0, 0.0, 1.0}, // blue
	{0.0, 1.0, 1.0}, // blue-sky
	{0.0, 1.0, 0.0}, // green
	{1.0, 1.0, 0.0}, // yellow
	{1.0, 0.0, 0.0} // red
};

class fColorPalette {
	int ndiv; ///< number of divisions
	int type; ///< Use PALETTE1, PALETTE2, ... ?
	std::vector<fColor> color_range;
public :
	fColorPalette(int _ndiv = 3, int _type = 1);
	fColor get_color(int i);
	int get_ncolors();
	std::vector<fColor> get_color_range();
};

#endif
