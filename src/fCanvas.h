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
	int label_size;
	int stick_size;
	int stick_width;
	int frame_line_width;
	
	bool draw_secondary_stick = true;	
	bool coord_system_not_defined = true;
	double linear_transformation(double x1, double y1, double x2, double y2, double x) const; ///< match [x1, x2] to [y1, y2] or ([y2, y1] if y2 < y1) f(x1) = y1 and f(x2) = y2, return y = f(x)
public :
	fCanvas(int width, int height, double xmin, double xmax, double ymin, double ymax);
	int x2w(double x) const; ///< convert x to width (pixel system)
	int y2h(double y) const; ///< convert y to height (pixel system)
	double w2x(double w) const; ///< convert width to x 
	double h2y(double h) const; ///< convert height to y
	void define_coord_system(const Cairo::RefPtr<Cairo::Context>& cr); ///< must be used only one time !
	void do_not_draw_secondary_stick();
	void draw_frame(const Cairo::RefPtr<Cairo::Context>& cr);
	void draw_title(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const;
	void draw_xtitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const;
	void draw_ytitle(const Cairo::RefPtr<Cairo::Context>& cr, std::string text) const;
	
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
	void set_label_size(double s);
	void set_stick_size(double s);
	void set_stick_width(double s);
	void set_frame_line_width(double s);

	int    get_top_margin() const;
	int    get_bottom_margin() const;
	int    get_left_margin() const;
	int    get_right_margin() const;
	int    get_weff() const;
	int    get_heff() const;
	int    get_seff() const;
	double get_x_start() const;
	double get_x_end() const;
	double get_y_start() const;
	double get_y_end() const;
	fAxis  get_x_axis() const;
	fAxis  get_y_axis() const;
	int    get_title_size() const;
	int    get_label_size() const;
	int    get_stick_size() const;
	int    get_stick_width() const;
	int    get_frame_line_width() const;

};


#endif
