/**************************************************
 * fCanvas.cpp
 *
 * Define a canvas in cairo context, set the axis
 * system.
 *
 * @author Felix Touchte Codjo
 * @date February 16, 2025
 * **********************************************/

#include "fCanvas.h"
#include <cstdlib>


fCanvas::fCanvas(int _width, int _height, double xmin, double xmax, double ymin, double ymax) {
	width = _width;
	height = _height;
	double plot_ratio = 0.10;
	int window_size = std::min(width, height);
	top_margin = plot_ratio*window_size; 
	bottom_margin = plot_ratio*window_size;
	left_margin = 1.5*plot_ratio*window_size;
	right_margin = 0.5*plot_ratio*window_size;
	weff = width - left_margin - right_margin;
	heff = height - top_margin - bottom_margin;
	seff = std::min(weff, heff);

	double margin_ratio = 0.05;	
	x_start = (xmin == 0) ? 0 : xmin - margin_ratio*(xmax - xmin);
	x_end = xmax + margin_ratio*(xmax - xmin);
	y_start = (ymin == 0) ? 0 : ymin - margin_ratio*(ymax - ymin);
	y_end = ymax + margin_ratio*(ymax - ymin);
	ax = fAxis(x_start, x_end, 10, 5, 0);
	ay = fAxis(y_start, y_end, 10, 5, 0);

	title_size = 0.4*top_margin;	
	label_size = 0.3*std::min(bottom_margin, left_margin);
        stick_size = 0.03*seff;
        stick_width = 0.005*seff;
        frame_line_width = 0.01*seff;
}

double fCanvas::linear_transformation(double x1, double y1, double x2, double y2, double x) const {
	if (x1 == x2) {return x;} // do nothing
	double slope = (y2-y1)/(x2-x1);
	double y=  slope*(x-x1) + y1;
	return y;
}

int fCanvas::x2w(double x) const {
	return linear_transformation(x_start, 0, x_end, weff, x);
}

int fCanvas::y2h(double y) const {
	return linear_transformation(y_start, 0, y_end, -heff, y); // minus heff because of the axis orientation
}

void fCanvas::define_coord_system(const Cairo::RefPtr<Cairo::Context>& cr) {
	if (coord_system_not_defined) {
		cr->translate(left_margin, top_margin + heff);
		coord_system_not_defined = false;
	}
}

void fCanvas::do_not_draw_secondary_stick(){
	draw_secondary_stick = false;
}

void fCanvas::draw_frame(const Cairo::RefPtr<Cairo::Context>& cr) {
	//cr->translate(left_margin, top_margin + heff);
	define_coord_system(cr);
	// Draw the frame for axis
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->set_line_width(frame_line_width);
	cr->rectangle(0,0,weff,-heff);
	cr->stroke();
	// Draw main sticks x
	for (std::string s : ax.get_labels1()) {
		double value = std::atof(s.c_str());
		if ((value >= x_start) && (value <= x_end)) {
			// draw stick
			cr->set_source_rgb(0.0, 0.0, 0.0);
			cr->set_line_width(stick_width);
			cr->move_to(x2w(value), 0);
			cr->line_to(x2w(value), -stick_size);
			cr->stroke();
			// draw label
			cr->set_source_rgb(0.0, 0.0, 0.0);
			cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
			cr->set_font_size(label_size);
			Cairo::TextExtents te;
        		cr->get_text_extents(s, te);
			cr->move_to(x2w(value) - 0.5*te.width, 0.1*bottom_margin + te.height);
			cr->show_text(s.c_str());
		}
	}
	// Draw main sticks y
	for (std::string s : ay.get_labels1()) {
		double value = std::atof(s.c_str());
		if ((value >= y_start) && (value <= y_end)) {
			// draw stick
			cr->set_source_rgb(0.0, 0.0, 0.0);
			cr->set_line_width(stick_width);
			cr->move_to(0, y2h(value));
			cr->line_to(stick_size, y2h(value));
			cr->stroke();
			// draw label
			cr->set_source_rgb(0.0, 0.0, 0.0);
			cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
			cr->set_font_size(label_size);
			Cairo::TextExtents te;
			cr->get_text_extents(s, te);
			cr->move_to(-left_margin*0.1 - te.width, y2h(value) + 0.5*te.height);
			cr->show_text(s.c_str());
		}
	}
	if (draw_secondary_stick) {
		// Draw secondary sticks x
		for (std::string s : ax.get_labels2()) {
			double value = std::atof(s.c_str());
			if ((value >= x_start) && (value <= x_end)) {
				//draw stick
				cr->set_source_rgb(0.0, 0.0, 0.0);
				cr->set_line_width(stick_width);
				cr->move_to(x2w(value), 0);
				cr->line_to(x2w(value), -0.7*stick_size);
				cr->stroke();
				// draw label
				/*cr->set_source_rgb(0.0, 0.0, 0.0);
				cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
				cr->set_font_size(label_size);
				Cairo::TextExtents te;
				cr->get_text_extents(s, te);
				//cr->move_to(0.5*weff - 0.5*te.width, -heff-0.2*top_margin);	
				cr->move_to(x2w(value) - 0.5*te.width, 0.1*bottom_margin + te.height);
				cr->show_text(s.c_str());*/
			}
		}
		// Draw seconday sticks y
		for (std::string s : ay.get_labels2()) {
			double value = std::atof(s.c_str());
			if ((value >= y_start) && (value <= y_end)) {
				// draw stick
				cr->set_source_rgb(0.0, 0.0, 0.0);
				cr->set_line_width(stick_width);
				cr->move_to(0, y2h(value));
				cr->line_to(0.7*stick_size, y2h(value));
				cr->stroke();
				// draw label
				/*cr->set_source_rgb(0.0, 0.0, 0.0);
				cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
				cr->set_font_size(label_size);
				Cairo::TextExtents te;
				cr->get_text_extents(s, te);
				cr->move_to(-left_margin*0.1 - te.width, y2h(value) + 0.5*te.height);
				cr->show_text(s.c_str());*/
			}
		}
	}
	///////////
}

void fCanvas::draw_title(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const {
	// draw label
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
	cr->set_font_size(title_size);
	Cairo::TextExtents te;
	cr->get_text_extents(text, te);
	cr->move_to(0.5*weff - 0.5*te.width, -heff-0.2*top_margin);
	cr->show_text(text.c_str());
}

void fCanvas::draw_xtitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const {
	// draw label
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
	cr->set_font_size(title_size);
	Cairo::TextExtents te;
	cr->get_text_extents(text, te);
	cr->move_to(weff - te.width, 0.95*bottom_margin);
	cr->show_text(text.c_str());
}

void fCanvas::draw_ytitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const {
	// draw label
	cr->save();
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
	cr->set_font_size(title_size);
	cr->rotate_degrees(-90);
	//cr->move_to(-left_margin*0.9, -heff);
	//x -> up and y -> right
	Cairo::TextExtents te;
	cr->get_text_extents(text, te);
	cr->move_to(heff - te.width, -0.97*left_margin + te.height);
	cr->show_text(text.c_str());
	cr->restore();
}

void fCanvas::set_top_margin(int margin) { 
	if (height - margin - bottom_margin < 0) {return ;}
	top_margin = margin;
	// update
	heff = height - top_margin - bottom_margin;
	seff = std::min(weff, heff);
	stick_size = 0.025*seff;
        stick_width = 0.005*seff;
        frame_line_width = 0.01*seff;
	title_size = 0.4*top_margin;
	label_size = 0.3*std::min(bottom_margin, left_margin);;
}

void fCanvas::set_bottom_margin(int margin) { 
	if (height - top_margin - margin < 0) {return ;}
	bottom_margin = margin;
	// update
	heff = height - top_margin - bottom_margin;
	seff = std::min(weff, heff);
	stick_size = 0.025*seff;
        stick_width = 0.005*seff;
        frame_line_width = 0.01*seff;
	label_size = 0.3*std::min(bottom_margin, left_margin);;
}

void fCanvas::set_left_margin(int margin) { 
	if (width - margin - right_margin < 0) {return ;}
	left_margin = margin;
	// update
	weff = width - left_margin - right_margin;
	seff = std::min(weff, heff);
	stick_size = 0.025*seff;
        stick_width = 0.005*seff;
        frame_line_width = 0.01*seff;
	label_size = 0.3*std::min(bottom_margin, left_margin);
}

void fCanvas::set_right_margin(int margin) { 
	if (width - left_margin - margin < 0) {return ;}
	right_margin = margin;
	// update
	weff = width - left_margin - right_margin;
	seff = std::min(weff, heff);
	stick_size = 0.025*seff;
        stick_width = 0.005*seff;
        frame_line_width = 0.01*seff;
}

void fCanvas::set_x_start(double value) { 
	x_start = value;
	// update
	ax = fAxis(x_start, x_end, 10, 0);
}

void fCanvas::set_x_end(double value) {
	x_end = value;
	// update
	ax = fAxis(x_start, x_end, 10, 0);
}

void fCanvas::set_y_start(double value) { 
	y_start = value;
	// update
	ay = fAxis(y_start, y_end, 10, 0);
}

void fCanvas::set_y_end(double value) { 
	y_end = value;
	// update
	ay = fAxis(y_start, y_end, 10, 0);
}

void fCanvas::set_x_axis(fAxis _ax) {
	ax = _ax;
	// update
	x_start = ax.get_start();
	y_end = ax.get_end();
}

void fCanvas::set_y_axis(fAxis _ay) { 
	ay = _ay;
	// update
	y_start = ay.get_start();
	y_end = ay.get_end();
}

void fCanvas::set_title_size(double s) { title_size = s*top_margin;}
void fCanvas::set_label_size(double s) { label_size = s*std::min(bottom_margin, left_margin);}
void fCanvas::set_stick_size(double s) { stick_size = s*seff;}
void fCanvas::set_stick_width(double s) { stick_width = s*seff;}
void fCanvas::set_frame_line_width(double s) { frame_line_width = seff*s;}

int    fCanvas::get_top_margin() const {return top_margin;}
int    fCanvas::get_bottom_margin() const {return bottom_margin;}
int    fCanvas::get_left_margin() const {return left_margin;}
int    fCanvas::get_right_margin() const {return right_margin;}
int    fCanvas::get_weff() const {return weff;}
int    fCanvas::get_heff() const {return heff;}
int    fCanvas::get_seff() const {return seff;}
double fCanvas::get_x_start() const {return x_start;}
double fCanvas::get_x_end() const {return x_end;}
double fCanvas::get_y_start() const {return y_start;}
double fCanvas::get_y_end() const {return y_end;}
fAxis  fCanvas::get_x_axis() const {return ax;}
fAxis  fCanvas::get_y_axis() const {return ay;}
int    fCanvas::get_title_size() const {return title_size;}
int    fCanvas::get_label_size() const {return label_size;}
int    fCanvas::get_stick_size() const {return stick_size;}
int    fCanvas::get_stick_width() const {return stick_size;}
int    fCanvas::get_frame_line_width() const {return frame_line_width;}


