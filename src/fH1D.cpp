/***********************************************
 * Class for 1D histogram
 *
 * designed to be used in gtkmm
 * drawing area
 *
 * @author Felix Touchte Codjo
 * @date February 12, 2025
 * ********************************************/
#include "fH1D.h"
#include "fCanvas.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <functional>

fH1D::fH1D(std::string _title, int _nbins, double _xmin, double _xmax) : title(_title), nbins(_nbins), xmin(_xmin), xmax(_xmax) {
	if (xmax < xmin) {
		printf("Histogram parameters are incorrects : xmax < xmin");		
		return ;
	}
	binw = (xmax - xmin)/nbins;
	for (int i = 0; i < nbins; i++) {
		binArray.push_back(xmin + i*binw + 0.5*binw);
		binBuffer.push_back(0.0);
	}
	underflow = 0;
	overflow = 0;
	nEntries = 0;
	sum = 0;
	sum2 = 0;
}

void fH1D::fill(double x) {
	int bin = getBinNumber(x);
	if (bin == -1) { underflow++; return;}
	if (bin == -11) { overflow++; return;}
	nEntries++;
	binBuffer.at(bin) = binBuffer.at(bin) + 1.0;
	// stats
	sumw += 1;
	sum += x;
	sum2 += x*x;
}

/**
 * fill x with weight of w
 * @param x value to fill
 * @param w weight
 *
 * @note return -1 if undeflow and -11 if overflow
 */
void fH1D::fill(double x, double w) {
	int bin = getBinNumber(x);
	if (bin == -1) { underflow++; return ;}
	if (bin == -11) { overflow++; return ;}
	nEntries++;
	binBuffer.at(bin) = binBuffer.at(bin) + w;
	//stats
	sumw += w;
	sum += w*x;
	sum2 += w*x*x;
}
/**
 *
 * @note numerotation starts at 0
 */
int fH1D::getBinNumber(double x) const {
	if (x < xmin) {return -1;}
	if (x > xmax) {return -11;}
	for (int i = 0; i < nbins; i++) {
		double xinf = xmin + i*binw;
		double xsup = xinf + binw;
		if ((x >= xinf) && (x < xsup)) {
			return i;
		}
	}
	return -1;
}

double fH1D::getBinBufferContent(int bin) const {
	if ((bin < 0) || (bin >= nbins)) {
		return 0;
	}
	return binBuffer.at(bin);
}

double fH1D::getBinArrayContent(int bin) const {
	/*if ((bin < 0) || (bin >= nbins)) {
		return 0;
	}*/
	return binArray.at(bin);
}

unsigned long int fH1D::getEntries() const { return nEntries;}
double fH1D::getMean() const { return sum/sumw;}
double fH1D::getStDev() const { return sqrt(sum2/sumw - getMean()*getMean());}
double fH1D::getBinWidth() const {return binw;}
int fH1D::getNumberOfBins() const {return nbins;}
std::vector<double> fH1D::getBinArray() const {return binArray;}
std::vector<double> fH1D::getBinBuffer() const {return binBuffer;}
void fH1D::set_xtitle(std::string name) {xtitle = name;}
void fH1D::set_ytitle(std::string name) {ytitle = name;}

double fH1D::getMax() const {
	int vmax = 0;
	for (int i = 0; i < nbins; i++) {
		vmax = (vmax < binBuffer.at(i)) ? binBuffer.at(i) : vmax;
	}
	return vmax;
}


void fH1D::reset() {
	if (binArray.size() < 1) {return ;}
	for (int i = 0; i < nbins; i++) {
                binBuffer.at(i) = 0.0;
        }
        underflow = 0;
        overflow = 0;
        nEntries = 0;
        sum = 0;
        sum2 = 0;	
}

void fH1D::print() {
	printf("Title : %s , nEntries : %ld , mean : %lf , stdev : %lf \n", title.c_str(), getEntries(), getMean(), getStDev());
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
}

//inclure tout gtkmm
//ou télécharger cairo
//draw(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height){
void fH1D::draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	// Define the main canvas
	fCanvas canvas(width, height, xmin, xmax, 0, getMax());
	canvas.define_coord_system(cr);
	//canvas.draw_frame(cr);
	canvas.draw_title(cr, title);
	canvas.draw_xtitle(cr, xtitle);
	canvas.draw_ytitle(cr, ytitle);
	// Draw contour	
	cr->set_source_rgb(0.0, 0.0, 1.0);
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
	cr->stroke_preserve(); // preserve the path
	cr->close_path();
	cr->set_source_rgb(fill_color.r, fill_color.g, fill_color.b);
	cr->fill();
	canvas.set_frame_line_width(0.005);
	canvas.draw_frame(cr); // draw frame and axis at the end

}

void fH1D::set_fill_color(fColor color) {
	fill_color = color;
}
