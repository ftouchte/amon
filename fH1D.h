/***********************************************
 * Class for 1D histogram
 *
 * designed to be used in gtkmm
 * drawing area
 *
 * @author Felix Touchte Codjo
 * @date February 12, 2025
 * ********************************************/

#ifndef F_HIST_1D
#define F_HIST_1D

#include <string>
#include <vector>
#include <gtkmm.h>

#include "fAxis.h"

struct fColor {
        double r;
        double g;
        double b;
};

class fH1D {
private :
	std::string title; ///< main annotation of the graph
	std::string xtitle; ///< name of the x axis
	std::string ytitle; ///< name of the y axis
	fColor fill_color = {1.0, 1.0, 1.0};
	int nbins; ///< number of bins	
	double xmin; ///< lower x value
	double xmax; ///< upper x value
	double binw; ///< bin width
	std::vector<double> binArray; ///< vector of centers of each bins
        std::vector<double> binBuffer; ///< occupancy of each bins
	// Statictics
	unsigned long int nEntries; ///< Number of entries without overflow and underflow
	int overflow; ///< number of overflow
	int underflow; ///< number of underflow
	double sumw; ///< sum of weight
	double sum; ///< sum of x without overflow and underflow
	double sum2; ///< sum of x*x without overflow and underflow


public :
	fH1D(std::string _title, int _nbins, double _xmin, double _xmax);
	void fill(double x);
	void fill(double x, double w);
	int getBinNumber(double x) const;
	double getBinArrayContent(int bin) const;
	double getBinBufferContent(int bin) const;
	unsigned long int getEntries() const;
	double getMean() const;
	double getStDev() const;
	double getBinWidth() const;
	int getNumberOfBins() const;
	std::vector<double> getBinArray() const;
	std::vector<double> getBinBuffer() const;
	void set_xtitle(std::string name);
	void set_ytitle(std::string name);
	double getMax() const;
	void reset();
	void draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void set_fill_color(fColor color);
	void print();
};


#endif
