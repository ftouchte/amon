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

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <functional>

fH2D::fH2D(std::string _title, int _nx, double _xmin, double _xmax, int _ny, double _ymin, double _ymax) 
	: title(_title), nx(_nx), xmin(_xmin), xmax(_xmax), ny(_ny), ymin(_ymin), ymax(_ymax), Hx("", _nx, _xmin, _xmax), Hy("", _ny, _ymin, _ymax)
{
	binwx = Hx.getBinWidth();
	binwy = Hy.getBinWidth();
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			binBuffer.push_back(0.0);
		}
	}
	overflow = 0;
	nEntries = 0;
}


void fH2D::fill(double x, double y) {
	Hx.fill(x);
	Hy.fill(y);
	int bin = getBinNumber(x,y);
	if (bin == -1) { overflow++; return;}
	nEntries++;
	binBuffer.at(bin) = binBuffer.at(bin) + 1.0;
}

/**
 * fill x,y with weight of w
 * @param x,y value to fill
 * @param w weight
 * @note return -1 if overflow
 */
void fH2D::fill(double x, double y, double w) {
	Hx.fill(x,w);
	Hy.fill(y,w);
	int bin = getBinNumber(x,y);
	if (bin == -1) { overflow++; return ;}
	nEntries++;
	binBuffer.at(bin) = binBuffer.at(bin) + w;
}

/**
 * @note numerotation starts at 0
 */
int fH2D::getBinNumber(double x, double y) const {
	if ((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax)) {return -1;}
	int i = Hx.getBinNumber(x);
	int j = Hy.getBinNumber(y);
	return j*nx + i;
}

double fH2D::getBinBufferContent(int bin) const {
	return binBuffer.at(bin);
}


unsigned long int fH2D::getEntries() const { return nEntries;}
std::vector<double> fH2D::getBinBuffer() const {return binBuffer;}
void fH2D::set_xtitle(std::string name) {xtitle = name;}
void fH2D::set_ytitle(std::string name) {ytitle = name;}

double fH2D::getMax() const {
	int vmax = 0;
	for (int i = 0; i < nx*ny; i++) {
		vmax = (vmax < binBuffer.at(i)) ? binBuffer.at(i) : vmax;
	}
	return vmax;
}


void fH2D::print() {
	printf("Title : %s , nEntries : %ld , mean_x : %lf , stdev_x : %lf, mean_y : %lf , stdev_y : %lf \n", title.c_str(), getEntries(), Hx.getMean(), Hx.getStDev(), Hy.getMean(), Hy.getStDev());
	/*
	fAxis ay(0, 100);
	printf("\033[32m");
	printf("           ");
	for (std::string s : ay.get_labels1()) {
		//int pos = std::atof(s.c_str());
		//int size = s.size();
		printf("%s ", s.c_str());
	}
	printf("  ====> this axis need to be fixed !!!!");
	printf("\n");
	for (int bin = 0; bin < nbins; bin++) {
		printf("\033[32m");
		int height = 100*getBinBufferContent(bin)/getMax();
		printf("%10.2lf ", getBinArrayContent(bin));
		//printf("%*.d ", (int) ceil(log10(nbins)), bin);
		for (int h = 0; h < height; h++) {
			printf(u8"█");
		}
		printf("\n");
	}
	printf("\033[37m");
	*/
}

void fH2D::draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	// Define the main canvas
	fCanvas canvas(width, height, xmin, xmax, ymin, ymax);
	canvas.draw_frame(cr);
	canvas.draw_title(cr, title);
	canvas.draw_xtitle(cr, xtitle);
	canvas.draw_ytitle(cr, ytitle);
	// Draw contour	
	/*cr->set_source_rgb(0.0, 0.0, 1.0);
	cr->set_line_width(0.008*canvas.get_seff());
	cr->move_to(canvas.x2w(xmin), canvas.y2h(0.0));
	for (int bin = 0; bin < nbins; bin++) {
		double x = binArray[bin] - 0.5*binw;
		double y = binBuffer[bin];
		//cr->move_to();
		cr->line_to(canvas.x2w(x), canvas.y2h(y));
		cr->line_to(canvas.x2w(x + binw), canvas.y2h(y));
	}
	cr->line_to(canvas.x2w(xmax), canvas.y2h(0.0));
	cr->stroke();*/
}


