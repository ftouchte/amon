/************************************************
 * fGtkmm is a class to customize some gtkmm's 
 * widgets
 *
 * main purpose : migrate from 4.0 to 3.0
 * (strange but it is a constraint)
 *
 * @author ftouchte
 * @date March 11, 2025
 * *********************************************/

#ifndef F_GTKMM_H
#define F_GTKMM_H


#include <gtkmm.h>
#include <vector>
#include <string>

#include "fH1D.h"
#include "AhdcDetector.h"

class fAreaH1D : public Gtk::DrawingArea {
	fH1D *ptrH1D;
public :
	fAreaH1D (fH1D* _ptr);
	~fAreaH1D();
	bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr);
};


class fAreaWaveform : public Gtk::DrawingArea {
	std::vector<double> vx, vy;
	std::string text;
	double leadingEdgeTime, timeOverThreshold, constantFractionTime, adcMax, adcOffset;
public :
	fAreaWaveform (const std::vector<double> & _vx, const std::vector<double> & _vy, const std::string & _text);
	void set_decoded_output(double _leadingEdgeTime, double _timeOverThreshold, double _constantFractionTime, double _adcMax, double _adcOffset);
	bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr);
	void cairo_plot_graph(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
};

class fAreaAhdcView : public Gtk::DrawingArea {
	AhdcDetector *ahdc;
	std::vector<AhdcWire> ListOfWires;
public :
	fAreaAhdcView();
	void set_ListOfWires(AhdcDetector * _ahdc, std::vector<AhdcWire> _ListOfWires);
	~fAreaAhdcView();
	bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr);
};



#endif
