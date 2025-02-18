/***********************************************
 * Class for 2D histogram
 *
 * designed to be used in gtkmm
 * drawing area || cairo context
 *
 * @author Felix Touchte Codjo
 * @date February 17, 2025
 * ********************************************/

#ifndef F_HIST_2D
#define F_HIST_2D

#include <string>
#include <vector>
#include <gtkmm.h>
#include "fH1D.h"

class fH2D {
private :
	std::string title; ///< main annotation of the graph
	std::string xtitle; ///< name of the x axis
	std::string ytitle; ///< name of the y axis
	std::vector<double> binBuffer; ///< occupancy of each bins
	int nx; ///< number of x bins
	double xmin; ///< lower x value
	double xmax; ///< upper x value
	int ny; ///< number of y bins
	double ymin; ///< lower y value
        double ymax; ///< upper y value
	double binwx; ///< x bin width
	double binwy; ///< y bin width
	// Statictics
	unsigned long int nEntries; ///< Number of entries without overflow and underflow
	int overflow; ///< number of overflow
	fH1D Hx; ///< histogram of x
	fH1D Hy; ///< histogram of y
public :
	fH2D(std::string _title, int _nx, double _xmin, double _xmax, int _ny, double _ymin, double _ymax);
	void fill(double x, double y);
	void fill(double x, double y, double w);
	int getBinNumber(double x, double y) const;
	double getBinBufferContent(int bin) const;
	fH1D getHx() const;
	fH1D getHy() const;
	double getMax() const;
	unsigned long int getEntries() const;
	std::vector<double> getBinBuffer() const;
	void set_xtitle(std::string name);
	void set_ytitle(std::string name);
	void draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void print();
};


#endif
