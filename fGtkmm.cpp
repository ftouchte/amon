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

#include "fGtkmm.h"
#include "fCanvas.h"

/** Constructor */
fAreaH1D::fAreaH1D (fH1D* _ptr) : Gtk::DrawingArea(), ptrH1D(_ptr) {}

/** Destructor */
fAreaH1D::~fAreaH1D(){
	if (!ptrH1D) { delete ptrH1D;}
}

/** override on_draw() */
bool fAreaH1D::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
	const Gtk::Allocation allocation = get_allocation();
	int width = allocation.get_width();
	int height = allocation.get_height();
	ptrH1D->draw_with_cairo(cr, width, height);
	return true;
}

/** Constructor */
fAreaWaveform::fAreaWaveform(const std::vector<double> & _vx, const std::vector<double> & _vy, const std::string & _text) : Gtk::DrawingArea(), vx(_vx), vy(_vy), text(_text) {}

void fAreaWaveform::set_decoded_output(double _leadingEdgeTime, double _timeOverThreshold, double _constantFractionTime, double _adcMax, double _adcOffset){
	leadingEdgeTime = _leadingEdgeTime;
	timeOverThreshold = _timeOverThreshold;
	constantFractionTime = _constantFractionTime;
	adcMax = _adcMax;
	adcOffset = _adcOffset;
}

/** override on_draw() */
bool fAreaWaveform::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
	const Gtk::Allocation allocation = get_allocation();
	int width = allocation.get_width();
	int height = allocation.get_height();
	this->cairo_plot_graph(cr, width, height);
	return true;
}

/** Draw waveform + decoded output */
void fAreaWaveform::cairo_plot_graph(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height){
	// Determine the min and the max of the data	
	int Npts = vx.size();
	if (Npts != (int) vy.size()) {return ;}
	double xmin = vx[0], xmax = vx[0];
	double ymin = vy[0], ymax = vy[0];
	for (int i = 0; i < Npts; i++){
		xmin = (xmin < vx[i]) ? xmin : vx[i];
 		xmax = (xmax > vx[i]) ? xmax : vx[i];
		ymin = (ymin < vy[i]) ? ymin : vy[i];
		ymax = (ymax > vy[i]) ? ymax : vy[i];
	}
	fCanvas canvas(width, height, xmin, xmax, ymin, ymax);
	canvas.define_coord_system(cr);
	canvas.draw_title(cr, "");
	canvas.draw_xtitle(cr, "bin");
	canvas.draw_ytitle(cr, "adc");
	// x coord to width
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	// y coord to height
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};

	int seff = canvas.get_seff();
	int heff = canvas.get_heff();
	int weff = canvas.get_weff();	
	// Draw points
	cr->set_source_rgb(0.0, 0.0, 1.0);
	cr->set_line_width(0.01*seff);
	cr->move_to(x2w(vx[0]),y2h(vy[0]));
	for (int i = 1; i < Npts; i++) {
		// draw a line between points i and i-1
		cr->line_to(x2w(vx[i]),y2h(vy[i]));
	}
	cr->stroke();
	
	// ___________________________
	// Show decoded values
	// ___________________________
	
	/*
	std::vector<short> samples;
	for (int i = 0; i < (int) vy.size(); i++) {
		samples.push_back((short) vy[i]);
	}
	decoder.adcOffset = (short) (samples[0] + samples[1] + samples[2] + samples[3] + samples[4])/5;
	std::map<std::string,double> output = decoder.extract(samples);
	//double timeMax = output["timeMax"];
	double leadingEdgeTime = output["leadingEdgeTime"];
	double timeOverThreshold = output["timeOverThreshold"];
	double constantFractionTime = output["constantFractionTime"];
	double adcMax = output["adcMax"];
	double adcOffset = output["adcOffset"];
	*/

	// Display leadingEdgeTime
	cr->set_source_rgb(0.0, 1.0, 0.0); // green
	cr->set_line_width(0.01*seff);
	cr->move_to(x2w(leadingEdgeTime),0);
	cr->line_to(x2w(leadingEdgeTime),-heff);
	cr->stroke();

	// Display constantFractionTime
	cr->set_source_rgb(1.0, 0.0, 0.0); // red
	cr->set_line_width(0.01*seff);
	cr->move_to(x2w(constantFractionTime),0);
	cr->line_to(x2w(constantFractionTime),-heff);
	cr->stroke();

	// Display timeOverThreshold
	cr->set_source_rgb(0.016, 0.925, 1); // bleu ciel
	cr->set_line_width(0.01*seff);
	cr->move_to(x2w(leadingEdgeTime), y2h(adcOffset + adcMax*0.5));
	cr->line_to(x2w(leadingEdgeTime + timeOverThreshold), y2h(adcOffset + adcMax*0.5));
	cr->line_to(x2w(leadingEdgeTime + timeOverThreshold), 0);
	//cr->line_to(x2w(leadingEdgeTime + timeOverThreshold), -heff);
	cr->stroke();

	// Display adcMax
	cr->set_source_rgb(1.0, 0.871, 0.016); // rose
	cr->set_line_width(0.01*seff);
	cr->move_to(0,y2h(adcOffset + adcMax));
	cr->line_to(weff, y2h(adcOffset + adcMax));
	//cr->line_to(x2w(timeMax), y2h(adcOffset + adcMax));
	//cr->line_to(x2w(timeMax), 0);
	cr->stroke();

	// Display the layer ID
	cr->set_source_rgb(1.0, 0.0, 0.0);
	cr->select_font_face("@cairo:sans-serif",Cairo::FontSlant::FONT_SLANT_NORMAL,Cairo::FontWeight::FONT_WEIGHT_NORMAL);
	cr->set_font_size(seff*0.1);
	cr->move_to(weff*0.7, -heff*0.8);
	cr->show_text(text.c_str());
	
	// draw frame and axis
	canvas.draw_frame(cr);	
}


/** Constructor */
fAreaAhdcView::fAreaAhdcView() : Gtk::DrawingArea() {}

/** Destructor */
fAreaAhdcView::~fAreaAhdcView() {
	if (!ahdc) { delete ahdc;}
}

/** */
void fAreaAhdcView::set_ListOfWires(AhdcDetector * _ahdc, std::vector<AhdcWire> _ListOfWires) {
	ahdc = _ahdc;
	ListOfWires = _ListOfWires;
}

/** Override on_draw() */
bool fAreaAhdcView::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
	const Gtk::Allocation allocation = get_allocation();
	int width = allocation.get_width();
	int height = allocation.get_height();
	fCanvas canvas(width, height, -80, 80, -80, 80);
	int window_size = std::min(width,height);
	canvas.set_top_margin(0.05*window_size);
	canvas.set_bottom_margin(0.05*window_size);
	canvas.set_left_margin(0.05*window_size);
	canvas.set_right_margin(0.05*window_size);
	canvas.set_frame_line_width(0.005);
	canvas.set_stick_width(0.002);
	canvas.set_label_size(0.5);
	canvas.set_title_size(0.6);
	canvas.set_x_start(-80);
	canvas.set_x_end(80);
	canvas.set_y_start(-80);
	canvas.set_y_end(80);
	canvas.define_coord_system(cr);
	canvas.do_not_draw_secondary_stick();
	canvas.draw_title(cr, "face connected to the electronics");
	canvas.draw_xtitle(cr, "");
	canvas.draw_ytitle(cr, "");
	canvas.draw_frame(cr);
	// x coord to width
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	// y coord to height
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};

	int seff = canvas.get_seff();
	int heff = canvas.get_heff();
	int weff = canvas.get_weff();
	double x_start = canvas.get_x_start();
	double x_end = canvas.get_x_end();
	double y_start = canvas.get_y_start();
	double y_end = canvas.get_y_end();
	
	// Draw AHDC geometry	
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire wire = *ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					// max radius == 68 cm
					// the distance between to wires of the last layer (radius == 68 cm) is
					// d = | exp((n+1)*theta) - exp(n*theta) | * radius
					// theta == 360°/99; 99 is the number of wires of this layer
					// so d = radius*sqrt( (cos(theta) -1)^2 + sin(theta)^2) 
					// so d = 4.315
					// marker radius should be < d/2, we take 2.0
					double marker_size = std::min(2.0*weff/(x_end-x_start), 2.0*heff/(y_end-y_start));
					cr->set_line_width(0.002*seff);
					cr->set_source_rgba(0.0, 0.0, 0.0,0.5);
					// if the previous reference point of the cairo context is not in the curve (when using cr->arc(...))
					// a straight line is added from the previous point to the place where the arc is started
					// solution : move_to the start place before setting the path
					cr->move_to(x2w(wire.top.x) + marker_size, y2h(wire.top.y)); // correspond to angle == 0 
					cr->arc(x2w(wire.top.x), y2h(wire.top.y) , marker_size, 0, 2*M_PI);
					cr->stroke();
					if ((wire.top.x == 0) && (wire.top.y < 0)) {
						// these wires hava a component id == 1 (it is the start of the numerotation)
						cr->set_line_width(0.002*seff);
						cr->set_source_rgba(0.0, 1.0, 0.0, 0.5);
						cr->move_to(x2w(wire.top.x) + marker_size, y2h(wire.top.y));
						cr->arc(x2w(wire.top.x), y2h(wire.top.y) , marker_size, 0, 2*M_PI);
						cr->stroke();
					}
				}
			}
		}
	}
	// Show activated wires
	for (AhdcWire wire : ListOfWires) {
		cr->set_source_rgb(1.0, 0.0, 0.0);
		double marker_size = std::min(2.0*weff/(x_end-x_start), 2.0*heff/(y_end-y_start));
		cr->arc(x2w(wire.top.x), y2h(wire.top.y) , marker_size, 0, 2*M_PI);
		cr->fill();
	}

	return true;
}
