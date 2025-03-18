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
#include "fH1D.h"

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
	Gtk::Grid Grid_waveformsPerLayer;
	Gtk::DrawingArea DrawingArea_event;
	Gtk::DrawingArea DrawingArea_test;
	// page 1
	Gtk::Box HBox_histograms;
	Gtk::Grid Grid_histograms;
	
	/******************
         * FOOTER
         * ***************/
	
	Gtk::Box HBox_footer;
	//Gtk::Box HBox_prev, HBox_next, HBox_pause, HBox_run, HBox_info,  HBox_hipo4, HBox_reset;
	Gtk::Button Button_settings, Button_prev, Button_next, Button_pause, Button_run, Button_hipo4, Button_reset;
	Gtk::Image img_settings, img_prev, img_next, img_pause, img_run, img_hipo4, img_reset;
	Gtk::Box HBox_info;
	Gtk::Label Label_info;
	Gtk::Label Label_header;
	
	/**********************
	 * Settings
	 * *******************/
	//Gtk::Window Window_settings;
	//Gtk::Box VBox_settings;
	//Gtk::Box HBox_Scale_adcMax, HBox_Scale_leadingEdgeTime_min, HBox_Scale_leadingEdgeTime_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_adcMax; 
	Glib::RefPtr<Gtk::Adjustment> Adjustment_leadingEdgeTime_min, Adjustment_leadingEdgeTime_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_timeOverThreshold_min, Adjustment_timeOverThreshold_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_timeMax_min, Adjustment_timeMax_max;
	Gtk::Scale Scale_adcMax;
	Gtk::Scale Scale_leadingEdgeTime_min, Scale_leadingEdgeTime_max;
	Gtk::Scale Scale_timeOverThreshold_min, Scale_timeOverThreshold_max;
	Gtk::Scale Scale_timeMax_min, Scale_timeMax_max;
	
	/**********************
	 *  DATA
	 * *******************/
	std::string filename;
	hipo::reader hipo_reader;
	hipo::banklist hipo_banklist;
	long unsigned int hipo_nEvent = 0;
	long unsigned int hipo_nEventMax = 1;
	int nWF; ///< number of waveforms	
	AhdcDetector *ahdc; ///< AHDC detector
	std::vector<AhdcWire> ListOfWires; ///< List of wires activated
	std::vector<std::string> ListOfWireNames; ///< List of wire names
	std::vector<std::vector<short>> ListOfSamples;  ///< List of samples
	std::vector<double> ListOfAdc;		///< List of adcMax
	std::vector< std::vector<std::vector<short>>  > ListOfSamplesPerLayer;
	AhdcExtractor decoder;
	double adcCut = 0;
	double leadingEdgeTime_min = 0;
	double leadingEdgeTime_max = 49;
	double timeOverThreshold_min = 0;
	double timeOverThreshold_max = 49;
	double timeMax_min = 0;
	double timeMax_max = 49;
	
	// Histograms
	bool is_paused = false;
	bool is_reset = false;
	fH1D hist1d_adcMax;
	fH1D hist1d_leadingEdgeTime;
	fH1D hist1d_timeOverThreshold;
	fH1D hist1d_timeMax;
	fH1D hist1d_adcOffset;
	fH1D hist1d_constantFractionTime;	

public :
	Window();
	~Window();

	bool dataEventAction();
	void drawWaveforms();
	void drawWaveformsPerLayer();
	void drawHistograms();
	
	// Signals
	void on_button_settings_clicked();
	void on_button_prev_clicked();
	void on_button_next_clicked();
	void on_button_pause_clicked();
	void on_button_run_clicked();
	void on_button_hipo4_clicked();
	void on_button_reset_clicked();
	void on_book_switch_page(Gtk::Widget * _page, guint page_num); // le template est imposé ! page_num == Pages.get_current_page()
	void on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog); // used in on_button_hipo4_clicked()
	void on_draw_event(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void on_draw_test(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void cairo_plot_graph(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height, std::vector<double> vx, std::vector<double> vy, std::string annotation);
};

#endif


