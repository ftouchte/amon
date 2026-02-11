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
#include <cmath>
#include <map>
#include "hipo4/reader.h"

#include "AhdcDetector.h"
#include "AhdcExtractor.h"
#include "fH1D.h"
#include "fH2D.h"

class Window : public Gtk::Window {
protected :
	// Widgets
	Gtk::Box VBox_main, VBox_header, VBox_body, VBox_footer;
	
	/******************  HEADER  ***************
	 * *****************************************/
	
	/******************  BODY  ***************
	 * *****************************************/
	
	Gtk::Notebook Book;
	Gtk::Box HBox_eventViewer;
	Gtk::Grid Grid_eventViewer, Grid_waveforms; // 2x1 and NxP 
	Gtk::Grid Grid_waveformsPerLayer;
	Gtk::DrawingArea DrawingArea_event;
	Gtk::DrawingArea DrawingArea_test;
	Gtk::Box HBox_histograms;
	Gtk::Grid Grid_histograms;
	Gtk::Grid Grid_occupancy;
	
	Gtk::Box HBox_3Dview;
	Gtk::Box VBox_3Dview_settings;
	Gtk::DrawingArea DrawingArea_3Dview;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_alpha; 
	Glib::RefPtr<Gtk::Adjustment> Adjustment_beta; 
	Glib::RefPtr<Gtk::Adjustment> Adjustment_gamma; 
	Gtk::Scale Scale_alpha;
	Gtk::Scale Scale_beta;
	Gtk::Scale Scale_gamma;
	double angle_alpha = 0.45*M_PI;
	double angle_beta = 0.1*M_PI;
	double angle_gamma = 0.1*M_PI;
	/******************  FOOTER  ***************
	 * *****************************************/
	
	Gtk::Box HBox_footer;
	//Gtk::Box HBox_prev, HBox_next, HBox_pause, HBox_run, HBox_info,  HBox_hipo4, HBox_reset;
	Gtk::Button Button_settings, Button_prev, Button_next, Button_pause, Button_run, Button_hipo4, Button_reset, Button_zpos, Button_restart_histo;
	Gtk::Image img_settings, img_prev, img_next, img_pause, img_run, img_hipo4, img_reset, img_zpos, img_restart_histo;
	Gtk::Box HBox_info;
	Gtk::Label Label_info;
	Gtk::Label Label_header;
	
	/******************  FOOTER/Settings button  ***************
	 * *********************************************************/
	
	// Decoding variables (must be changed if changed in caotjava)	
	const int ADC_LIMIT;
	const int NumberOfBins;
	const double samplingTime;
	const double amplitudeFractionCFA;
	const int binDelayCFD;
	const double fractionCFD;

	Glib::RefPtr<Gtk::Adjustment> Adjustment_amplitude_min, Adjustment_amplitude_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_cut_adcOffset_min, Adjustment_cut_adcOffset_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_cut_leadingEdgeTime_min, Adjustment_cut_leadingEdgeTime_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_cut_timeOverThreshold_min, Adjustment_cut_timeOverThreshold_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_cut_timeMax_min, Adjustment_cut_timeMax_max;
	Glib::RefPtr<Gtk::Adjustment> Adjustment_zpos; 
	Gtk::Scale Scale_amplitude_min, Scale_amplitude_max;
	Gtk::Scale Scale_cut_adcOffset_min, Scale_cut_adcOffset_max;
	Gtk::Scale Scale_cut_leadingEdgeTime_min, Scale_cut_leadingEdgeTime_max;
	Gtk::Scale Scale_cut_timeOverThreshold_min, Scale_cut_timeOverThreshold_max;
	Gtk::Scale Scale_cut_timeMax_min, Scale_cut_timeMax_max;
	Gtk::Scale Scale_zpos;
	Gtk::CheckButton CheckButton_apply_wfType_cuts; 
	Gtk::CheckButton CheckButton_apply_raw_cuts; 
	Gtk::CheckButton CheckButton_wfType_0; 
	Gtk::CheckButton CheckButton_wfType_1; 
	Gtk::CheckButton CheckButton_wfType_2; 
	Gtk::CheckButton CheckButton_wfType_3; 
	Gtk::CheckButton CheckButton_wfType_4; 
	Gtk::CheckButton CheckButton_wfType_5; 
		
	/******************  DATA  *****************
	 * *****************************************/
	
	std::string filename;
	std::string pdf_output_name;
	hipo::reader hipo_reader;
	hipo::dictionary hipo_factory;
	hipo::event hipo_event;
	hipo::bank adcBank;
	hipo::bank wfBank;
	hipo::bank trackBank;
	long unsigned int hipo_nEvent = 0;
	long unsigned int hipo_nEventMax = 1;
	int nWF; ///< number of waveforms	
	int ntracks; ///< number of tracks	
	AhdcDetector *ahdc; ///< AHDC detector
	std::vector<double> ListOfAdc;		///< List of adcMax
    int minADC = 0;
    int maxADC = 4095;
	// Define variable for cut
	double cut_amplitude_min = 0;
	double cut_amplitude_max = 4095;
	double cut_adcOffset_min = 0;
	double cut_adcOffset_max = 1000;
	double cut_leadingEdgeTime_min = 0;
	double cut_leadingEdgeTime_max = 950;
	double cut_timeOverThreshold_min = 0;
	double cut_timeOverThreshold_max = 950;
	double cut_timeMax_min = 0;
	double cut_timeMax_max = 950;
	bool flag_apply_wfType_cuts = false;
	bool flag_apply_raw_cuts = false;
	bool flag_wfType_0 = false;
	bool flag_wfType_1 = false;
	bool flag_wfType_2 = false;
	bool flag_wfType_3 = false;
	bool flag_wfType_4 = false;
	bool flag_wfType_5 = false;
	double zpos = -150.0; ///< define in which z we should look at for AHDC (x,y) view
	int HV_SECTOR = -1;
	Gtk::TextView TextView_occupancy;
	Glib::RefPtr<Gtk::TextBuffer> TextBuffer_occupancy;
	// Histograms
	bool is_paused = false;
	bool is_reset = false;
	fH1D hist1d_adcMax;
	fH1D hist1d_leadingEdgeTime;
	fH1D hist1d_timeOverThreshold;
	fH1D hist1d_timeMax;
	fH1D hist1d_adcOffset;
	fH1D hist1d_constantFractionTime;
	fH2D hist2d_occupancy;
	std::map<std::string, std::function<void(const Cairo::RefPtr<Cairo::Context>&, int, int)>> renderers;

public :
	Window();
	~Window();

	bool dataEventAction();
	void drawWaveforms();
	void drawWaveformsPerLayer();
	void drawHistograms();
	void drawOccupancy();
	AhdcWire* getNearestAhdcWire(double x, double y, int & _layer, int & _component);
	void clearAhdcData();
	int getNumberOfWaveforms();
	bool is_oscillating(std::vector<double> samples);
	void getStats(double & MIN_ADC, double & MAX_ADC, int & MIN_OCC, int & MAX_OCC);
    bool update_cut_flag(double adcMax, double adcOffset, double leadingEdgeTime, double timeOverThreshold, double timeMax, int wfType);
    void update_gui();
    void restart_histograms();
    void updateWireStatus();
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
	void on_draw_occupancy(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void on_draw_test(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	void on_mouse_clicked(int n_press, double x, double y);
	void ask_user_confirmation(std::string name, int width, int height);
	void select_file_to_save_as_pdf(std::string name, int width, int height);
	void export_as_pdf(std::string name, int width, int height);
	//void on_question_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result);
	void on_zpos_value_changed();
	void on_button_zpos_clicked();
	void cairo_plot_waveform(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height, AhdcWire* wire, std::string annotation);
	void on_draw_3Dview(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	int layer2number(int digit);
	void Get_HV_sector(int sector, int layer, int component, int & crate, int & slot, int & channel, int & hv, int & sub_hv);
};

#endif


