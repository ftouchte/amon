/**************************************
 * GUI for ALERT monitoring
 * 
 * @author ftouchte
 * @date Jan 14, 2025
 * **********************************/

#ifndef AHDC_GUI
#define AHDC_GUI

#include <gtkmm.h>
#include <string>
#include <vector>
#include "reader.h"

#include "AhdcDetector.h"
#include "AhdcExtractor.h"

class Window : public Gtk::Window {
protected :
	// Widgets
	Gtk::Box VBox_main, VBox_header, VBox_body, VBox_footer;
	
	/******************
         * HEADER
         * ***************/


	/******************
         * BODY
         * ***************/
	
	Gtk::Notebook Book;
	// page 0
	Gtk::Box HBox_eventViewer;
	Gtk::Grid Grid_eventViewer, Grid_waveforms; // 2x1 and NxP 
	Gtk::Picture Picture_event;
	Gtk::DrawingArea DrawingArea_event;
	Gtk::DrawingArea DrawingArea_test;
	// page 1
	Gtk::Box HBox_histograms;

	/******************
         * FOOTER
         * ***************/
	Gtk::Box HBox_footer;
	//Gtk::Box HBox_prev, HBox_next, HBox_pause, HBox_run, HBox_info,  HBox_hipo4, HBox_reset;
	Gtk::Button Button_prev, Button_next, Button_pause, Button_run, Button_hipo4, Button_reset;
	Gtk::Image img_prev, img_next, img_pause, img_run, img_hipo4, img_reset;
	Gtk::Box HBox_info;
	Gtk::Label Label_info;

	/// DATA
	std::string filename;
	hipo::reader hipo_reader;
	hipo::banklist hipo_banklist;
	long unsigned int hipo_nEvent = 0;
	int nWF; ///< number of waveforms	
	AhdcDetector *ahdc; ///< AHDC detector
	std::vector<AhdcWire> ListOfWires; ///< List of wires activated
	std::vector<std::string> ListOfWireNames; ///< List of wire names
	std::vector<std::vector<short>> ListOfSamples;  ///< List of samples
	AhdcExtractor decoder;

	// Histograms
	

public :
	Window();
	~Window();

	void dataEventAction();
	
	// Signals
	void on_button_prev_clicked();
	void on_button_next_clicked();
	void on_button_pause_clicked();
	void on_button_run_clicked();
	void on_button_hipo4_clicked();
	void on_button_reset_clicked();
	void on_book_switch_page(Gtk::Widget * _pages, guint page_num); // le template est imposé ! page_num == Pages.get_current_page()
	void on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog); // used in on_button_hipo4_clicked()
	void on_draw_event(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void on_draw_test(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void cairo_plot_graph(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height, std::vector<double> vx, std::vector<double> vy, std::string annotation);
};

#endif


