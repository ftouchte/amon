/**************************************
 * f3Dutils.h
 *
 * @author Felix Touchte Codjo
 * @date April 6, 2025
 * ***********************************/

#include "f3Dutils.h"
#include "fCanvas.h"
#include <cmath>

Frame::Frame(double _alpha, double _beta, double _gamma) : alpha(_alpha), beta(_beta), gamma(_gamma) {
	update();
}

PseudoV3D Frame::getXAxis() { return vx;}
PseudoV3D Frame::getYAxis() { return vy;}
PseudoV3D Frame::getZAxis() { return vz;}

PseudoV3D Frame::get_PseudoV3D(double x, double y, double z) {
	double X = x*vx.X + y*vy.X + z*vz.X;
	double Y = x*vx.Y + y*vy.Y + z*vz.Y;
	return PseudoV3D({X,Y});
}
void Frame::update() {
	// vx
	vx.X = cos(alpha);
	vx.Y = sin(alpha);
	// vy
	vy.X = cos(alpha + beta);
	vy.Y = sin(alpha + beta);
	// vz
	vz.X = cos(gamma);
	vz.Y = sin(gamma);
}
void Frame::set_alpha(double _alpha) { alpha = _alpha; update(); }
void Frame::set_beta (double _beta)  { beta = _beta; update(); }
void Frame::set_gamma(double _gamma) { gamma = _gamma; update(); }

void Frame::test1(double alpha, double beta, double gamma, const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	fCanvas canvas(width, height, -150, 150, -150, 150);
	int window_size = std::min(width,height);
	//canvas.set_top_margin(0.05*window_size);
	//canvas.set_bottom_margin(0.20*window_size);
	//canvas.set_left_margin(0.05*window_size);
	//canvas.set_right_margin(0.15*window_size);
	canvas.set_frame_line_width(0.005);
	//canvas.set_stick_width(0.002);
	//canvas.set_label_size(0.5);
	//canvas.set_title_size(0.6);
	//canvas.set_x_start(-80);
	//canvas.set_x_end(80);
	//canvas.set_y_start(-80);
	//canvas.set_y_end(80);
	canvas.define_coord_system(cr);
	canvas.do_not_draw_secondary_stick();
	canvas.draw_title(cr, "");
	canvas.draw_xtitle(cr, "");
	canvas.draw_ytitle(cr, "");
	
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
	int marker_size = 4;

	/*************************************************************
	 **********************  Simu 3D view  **********************/
	
	// Draw center : stationary point
	/*cr->set_source_rgb(0.0, 0.0, 1.0);
	cr->move_to(x2w(0) + marker_size, y2h(0));
	cr->arc(x2w(0), y2h(0), marker_size, 0, 2*M_PI); 
	cr->fill();*/
		
	// Draw stationary x axis
	cr->set_line_width(0.001*seff); 
	cr->set_source_rgba(0.0, 0.0, 0.0, 0.1);
	cr->move_to(x2w(0), 0);
	cr->line_to(x2w(0), -heff); 
	cr->stroke();

	// Draw stationary y axis 
	cr->set_line_width(0.001*seff);
	cr->set_source_rgba(0.0, 0.0, 0.0, 0.1);
	cr->move_to(0, y2h(0));
	cr->line_to(weff, y2h(0)); 
	cr->stroke();

	// Frame
	Frame frame(alpha, beta, gamma);
	//Frame frame(0.05*M_PI, 0.1*M_PI, 0.5*M_PI);
	
	// Draw x axis
	PseudoV3D ux = frame.getXAxis();	
	cr->set_line_width(0.005*seff);
	cr->set_source_rgba(0.0, 0.0, 0.0, 1.0);
	cr->move_to(x2w(-130*ux.X), y2h(-130*ux.Y));
	cr->line_to(x2w(130*ux.X), y2h(130*ux.Y)); 
	cr->stroke();

	// Draw y axis
	PseudoV3D uy = frame.getYAxis();	
	cr->set_line_width(0.005*seff);
	cr->set_source_rgba(0.0, 0.0, 0.0, 1.0);
	cr->move_to(x2w(-130*uy.X), y2h(-130*uy.Y));
	cr->line_to(x2w(130*uy.X), y2h(130*uy.Y)); 
	cr->stroke();

	// Draw z axis
	PseudoV3D uz = frame.getZAxis();	
	cr->set_line_width(0.005*seff);
	cr->set_source_rgba(0.0, 0.0, 0.0, 1.0);
	cr->move_to(x2w(0.0*uz.X), y2h(0.0*uz.Y));
	cr->line_to(x2w(130*uz.X), y2h(130*uz.Y)); 
	cr->stroke();

	// Draw a circle in the (x,y) plan
	for (double z : std::vector<double>({0.0, 100.0})) {
		cr->set_line_width(0.002*seff);
		cr->set_source_rgba(1.0, 0.0, 0.0, 1.0);
		int Npts = 50;
		double radius = 50;
		for (int i = 0; i <= Npts; i++) {
			double angle = i*2*M_PI/Npts;
			double x = radius*cos(angle);
			double y = radius*sin(angle);
			PseudoV3D pos = frame.get_PseudoV3D(x, y, z);
			if (i == 0) {
				cr->move_to(x2w(pos.X), y2h(pos.Y));
			}
			else {
				cr->line_to(x2w(pos.X), y2h(pos.Y));
			}
		}
		cr->stroke();
	}
	{ // Draw an helix
		cr->set_line_width(0.002*seff);
		cr->set_source_rgba(0.0, 1.0, 0.0, 1.0);
		int Npts = 150;
		int Ncycle = 10;
		double radius = 50;
		for (int i = 0; i <= Npts; i++) {
			double angle = i*Ncycle*M_PI/Npts;
			double x = radius*cos(angle);
			double y = radius*sin(angle);
			PseudoV3D pos = frame.get_PseudoV3D(x, y, 100.0*i/Npts);
			if (i == 0) {
				cr->move_to(x2w(pos.X), y2h(pos.Y));
			}
			else {
				cr->line_to(x2w(pos.X), y2h(pos.Y));
			}
		}
		cr->stroke();
	}
}
