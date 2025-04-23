/***********************************************
 * Class for 2D histogram
 *
 * designed to be used in gtkmm
 * drawing area || cairo context
 *
 * @author Felix Touchte Codjo
 * @date February 12, 2025
 * ********************************************/

#include "fH2D.h"
#include "fAxis.h"
#include "fCanvas.h"
#include "fColorPalette.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <functional>

fH2D::fH2D(std::string _title, int _nx, double _xmin, double _xmax, int _ny, double _ymin, double _ymax) 
	: title(_title), nx(_nx), xmin(_xmin), xmax(_xmax), ny(_ny), ymin(_ymin), ymax(_ymax) 
{
	  thisto = new TH2D(title.c_str(), title.c_str(), nx, xmin, xmax, ny, ymin, ymax);
}

fH2D::~fH2D() {
	delete thisto;
}

void fH2D::reset() {
	thisto->Reset();
}
void fH2D::fill(double x, double y) {
	thisto->Fill(x,y);
}

void fH2D::fill(double x, double y, double w) {
	thisto->Fill(x,y,w);
}

TH2D* fH2D::get_root_histo() {
	return thisto;
}

void fH2D::set_xtitle(std::string name) {xtitle = name;}
void fH2D::set_ytitle(std::string name) {ytitle = name;}

void fH2D::draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	// Define the main canvas
	fCanvas canvas(width, height, xmin, xmax, ymin, ymax);
	canvas.set_x_start(xmin);
	canvas.set_x_end(xmax);
	canvas.set_y_start(ymin);
	canvas.set_y_end(ymax);
	canvas.define_coord_system(cr);
	//canvas.draw_frame(cr);
	canvas.draw_title(cr, title);
	canvas.draw_xtitle(cr, xtitle);
	canvas.draw_ytitle(cr, ytitle);
	// x coord to width
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	// y coord to height
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};
	// Find max value
	double vmax = 0.0;
	for (int ix = 1; ix <= nx; ix++) {
		for (int iy = 1; iy <= ny; iy++) {
			double value = thisto->GetBinContent(ix, iy);
			if (value > vmax) vmax = value;
		}
	}
	// Draw each (x,y) bins
	fColorPalette Palette(5, 2);
	int ncolors = Palette.get_ncolors();
	double wx = (xmax-xmin)/nx;
	double wy = (ymax-ymin)/ny;
	for (int ix = 1; ix <= nx; ix++) {
		for (int iy = 1; iy <= ny; iy++) {
			double value = thisto->GetBinContent(ix, iy); // here, the numeratation starts at 1
			if (value > 0) {
				double slope = (1.0*(ncolors-1))/(vmax - 0.0);
				int color = floor(slope*(value-0.0) + 0.0);
				//printf("value : %lfi, color : %d\n", value, color);
				fColor rgb = Palette.get_color(color);
				cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
				double x = xmin + ix*wx;
				double y = ymin + iy*wy;
				cr->move_to(x2w(x-wx), y2h(y-wy));
				cr->line_to(x2w(x), y2h(y-wy));
				cr->line_to(x2w(x), y2h(y));
				cr->line_to(x2w(x-wx), y2h(y));
				cr->close_path();
				cr->fill();
			}
		}
	}
	// Draw frame
	canvas.set_frame_line_width(0.005);
	canvas.draw_frame(cr); // draw frame and axis at the end
}


