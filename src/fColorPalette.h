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


// kBird from CERN ROOT Project
// https://github.com/root-project/root/blob/9bb2f5076762c4171bfa71f44da3ee74fbd27b5b/core/base/src/TColor.cxx#L3055
inline const std::vector<fColor> PALETTE2 = {
	// {r, g, b}
	{0.2082, 0.1664, 0.5293},
	{0.0592, 0.3599, 0.8684},
	{0.0780, 0.5041, 0.8385},
	{0.0232, 0.6419, 0.7914},
	{0.1802, 0.7178, 0.6425},
	{0.5301, 0.7492, 0.4662},
	{0.8186, 0.7328, 0.3499},
	{0.9956, 0.7862, 0.1968},
	{0.9764, 0.9832, 0.0539}
};

inline const std::vector<fColor> PALETTE3 = {
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
