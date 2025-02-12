/***********************************************
 * Class for 1D histogram
 *
 * designed to be used in gtkmm
 * drawing area
 *
 * @author Felix Touchte Codjo
 * ********************************************/

#ifndef F_HIST_1D
#define F_HIST_1D

#include <string>
#include <vector>
#include "fAxis.h"


class fH1D {
private :
	std::string title; ///< main annotation of the graph
	std::string xtitle; ///< name of the x axis
	std::string ytitle; ///< name of the y axis
	fAxis ax; ///< x axis
	fAxis ay; ///< y axis
	std::vector<double> binArray; ///< vector of centers of each bins
	std::vector<double> binBuffer; ///< occupancy of each bins
	int nbins; ///< number of bins	
	double xmin; ///< lower x value
	double xmax; ///< upper x value
	double binw; ///< bin width
	// Statictics
	unsigned long int nEntries; ///< Number of entries without overflow and underflow
	int overflow; ///< number of overflow
	int underflow; ///< number of underflow
	double sumw; ///< sum of weight
	double sum; ///< sum of x without overflow and underflow
	double sum2; ///< sum of x*x without overflow and underflow


public :
	fHist1D(std::string _title = "", int _nbins, double _xmin, double _xmax);
	void fill(double x);
	void fill(double x, double w);
	int getBinNumber(double x) const;
	double getBinValue(int bin) const;
	double getBinBuffer(int bin) const;
	int getEntries() const;
	double getMean() const;
	double getStDev() const;
	double getBinWidth() const;
	int getNumberOfBins() const;
	std::vector<double> getBinArray() const;
	std::vector<double> getBinBuffer() const;
	void set_xtitle(std::string name);
	void set_ytitle(std::string name);
	double getMax() const;
	void draw();
};


#endif
