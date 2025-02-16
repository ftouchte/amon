/**************************************************
 * fCanvas.h
 *
 * Define a canvas in cairo context, set the axis
 * system.
 *
 * @author Felix Touchte Codjo
 * @date February 16, 2025
 * **********************************************/

#ifndef F_CANVAS_H
#define F_CANVAS_H

#include "fAxis.h"
#include <gtkmm.h>
#include <string>

class fCanvas {
private :
	int width;
	int height;
	int top_margin;
	int bottom_margin;
	int left_margin;
	int right_margin;
	int weff; ///< effective width size
	int heff; ///< effective height size
	int seff; ///< min of weff and heff
	
	double x_start;
	double x_end;
	double y_start;
	double y_end;
	fAxis ax; ///< X axis  
	fAxis ay; ///< Y axis

	int title_size;	
	int xlabel_size;
	int ylabel_size;
	int stick_size;
	int stick_width;
	int frame_line_width;

	double linear_transformation(double x1, double y1, double x2, double y2, double x); ///< match [x1, x2] to [y1, y2] or ([y2, y1] if y2 < y1) f(x1) = y1 and f(x2) = y2, return y = f(x)
public :
	fCanvas(int width, int height, double xmin, double xmax, double ymin, double ymax);
	int x2w(double x); ///< convert x to width (pixel system)
	int y2h(double y); ///< convert y to height (pixel system)
	//void update_axis_limit(); ///< update x_start, x_end, y_start, y_end
	void draw_frame(const Cairo::RefPtr<Cairo::Context>& cr);
	void draw_title(const Cairo::RefPtr<Cairo::Context>& cr, std::string text);
	void draw_xtitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text);
	void draw_ytitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text);
	
	void set_top_margin(int margin);
	void set_bottom_margin(int margin);
	void set_left_margin(int margin);
	void set_right_margin(int margin);
	void set_x_start(double value);
	void set_x_end(double value);
	void set_y_start(double value);
	void set_y_end(double value);
	void set_x_axis(fAxis _ax);
	void set_y_axis(fAxis _ay);
	void set_title_size(double s);
	void set_xlabel_size(double s);
	void set_ylabel_size(double s);
	void set_stick_size(double s);
	void set_stick_width(double s);
	void set_frame_line_width(double s);

	int    get_top_margin();
	int    get_bottom_margin();
	int    get_left_margin();
	int    get_right_margin();
	int    get_weff();
	int    get_heff();
	double get_x_start();
	double get_x_end();
	double get_y_start();
	double get_y_end();
	fAxis  get_x_axis();
	fAxis  get_y_axis();
	int    get_title_size();
	int    get_xlabel_size();
	int    get_ylabel_size();
	int    get_stick_size();
	int    get_stick_width();
	int    get_frame_line_width();

};


#endif
