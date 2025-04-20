/***********************************************
 * Class for 2D histogram
 *
 * designed to be used in gtkmm
 * drawing area || cairo context
 * 
 * wrapper of the ROOT TH2D
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
#include "TH2D.h"

class fH2D {
private :
	std::string title; ///< main annotation of the graph
	std::string xtitle; ///< name of the x axis
	std::string ytitle; ///< name of the y axis
	int nx; ///< number of x bins
	double xmin; ///< lower x value
	double xmax; ///< upper x value
	int ny; ///< number of y bins
	double ymin; ///< lower y value
        double ymax; ///< upper y value
	TH2D* thisto = nullptr; ///< ROOT 2D histograms
public :
	fH2D(std::string _title, int _nx, double _xmin, double _xmax, int _ny, double _ymin, double _ymax);
	~fH2D();
	void reset();
	void fill(double x, double y);
	void fill(double x, double y, double w);
	TH2D* get_root_histo();
	void set_xtitle(std::string name);
	void set_ytitle(std::string name);
	void draw_with_cairo(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
};


#endif
