/**************************************
 * GUI for ALERT monitoring
 * 
 * @author ftouchte
 * @date Jan 14, 2025
 * **********************************/

#include "gui.h"
#include "fAxis.h"
#include "fCanvas.h"
#include "fColorPalette.h"
#include "AhdcMapping.h"
#include "f3Dutils.h"

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <functional>
#include <cstdio>
#include <cstdlib>

#include "TString.h"
#include "TCanvas.h"

/** Constructor */
Window::Window() :
	// Initialisation // Take of the order of declaration
	VBox_main(Gtk::Orientation::VERTICAL),
	VBox_header(Gtk::Orientation::VERTICAL,10),
	VBox_body(Gtk::Orientation::VERTICAL,10),
	VBox_footer(Gtk::Orientation::VERTICAL,10),
	HBox_eventViewer(Gtk::Orientation::HORIZONTAL,10),
	HBox_histograms(Gtk::Orientation::HORIZONTAL,10),
	HBox_3Dview(Gtk::Orientation::HORIZONTAL,10),
	VBox_3Dview_settings(Gtk::Orientation::VERTICAL,10),
	// 3D view settings
	Adjustment_alpha(Gtk::Adjustment::create(0.45*M_PI, 0.0, 2*M_PI, 0.01, 0.0, 0.0)),
	Adjustment_beta(Gtk::Adjustment::create(0.1*M_PI, 0.0, 2*M_PI, 0.01, 0.0, 0.0)),
	Adjustment_gamma(Gtk::Adjustment::create(0.1*M_PI, 0.0, 2*M_PI, 0.01, 0.0, 0.0)),
	Scale_alpha(Adjustment_alpha, Gtk::Orientation::HORIZONTAL),
	Scale_beta(Adjustment_beta, Gtk::Orientation::HORIZONTAL),
	Scale_gamma(Adjustment_gamma, Gtk::Orientation::HORIZONTAL),
	// Footer
	HBox_footer(Gtk::Orientation::HORIZONTAL,10),
	HBox_info(Gtk::Orientation::HORIZONTAL,15),
	// Decoding parameters	
	ADC_LIMIT(4095),
	NumberOfBins(20),
	samplingTime(50.0), // set at 50.0 if you want use measure in terms of bins
	amplitudeFractionCFA(0.5),
	binDelayCFD(5),
	fractionCFD(0.3),
	// Value, lower, upper, step_increment, page_increment, page_size
	Adjustment_adcMax(Gtk::Adjustment::create(0.0, 0.0, ADC_LIMIT, 1, 0.0, 0.0)),
	Adjustment_cut_adcOffset_min(Gtk::Adjustment::create(0.0, 0.0, 1000, 1, 0.0, 0.0)),
	Adjustment_cut_adcOffset_max(Gtk::Adjustment::create(ADC_LIMIT, 0.0, 1000, 1, 0.0, 0.0)),
	Adjustment_cut_leadingEdgeTime_min(Gtk::Adjustment::create(0.0, 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_cut_leadingEdgeTime_max(Gtk::Adjustment::create(samplingTime*(NumberOfBins-1), 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_cut_timeOverThreshold_min(Gtk::Adjustment::create(0.0, 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_cut_timeOverThreshold_max(Gtk::Adjustment::create(samplingTime*(NumberOfBins-1), 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_cut_timeMax_min(Gtk::Adjustment::create(0.0, 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_cut_timeMax_max(Gtk::Adjustment::create(samplingTime*(NumberOfBins-1), 0.0, samplingTime*(NumberOfBins-1), 1, 0.0, 0.0)),
	Adjustment_zpos(Gtk::Adjustment::create(-150, -150.0, 150.0, 1.0, 0.0, 0.0)),
	Scale_adcMax(Adjustment_adcMax, Gtk::Orientation::HORIZONTAL),
	Scale_cut_adcOffset_min(Adjustment_cut_adcOffset_min, Gtk::Orientation::HORIZONTAL),
	Scale_cut_adcOffset_max(Adjustment_cut_adcOffset_max, Gtk::Orientation::HORIZONTAL),
	Scale_cut_leadingEdgeTime_min(Adjustment_cut_leadingEdgeTime_min, Gtk::Orientation::HORIZONTAL),
	Scale_cut_leadingEdgeTime_max(Adjustment_cut_leadingEdgeTime_max, Gtk::Orientation::HORIZONTAL),
	Scale_cut_timeOverThreshold_min(Adjustment_cut_timeOverThreshold_min, Gtk::Orientation::HORIZONTAL),
	Scale_cut_timeOverThreshold_max(Adjustment_cut_timeOverThreshold_max, Gtk::Orientation::HORIZONTAL),
	Scale_cut_timeMax_min(Adjustment_cut_timeMax_min, Gtk::Orientation::HORIZONTAL),
	Scale_cut_timeMax_max(Adjustment_cut_timeMax_max, Gtk::Orientation::HORIZONTAL),
	Scale_zpos(Adjustment_zpos, Gtk::Orientation::HORIZONTAL),
	// histograms
	hist1d_adcMax("adcMax", 200, 0.0, 4095.0),
        hist1d_leadingEdgeTime("leadingEdgeTime", 100, 0.0, samplingTime*(NumberOfBins-1)),
        hist1d_timeOverThreshold("timeOverThreshold", 100, 0.0, samplingTime*(NumberOfBins-1)),
        hist1d_timeMax("timeMax", 100, 0.0, samplingTime*(NumberOfBins-1)),
	hist1d_adcOffset("adcOffset (s1+s2+s3+s4+s5)/5", 100, 0.0, 1000),
        hist1d_constantFractionTime("constantFractionTime", 100, 0.0, samplingTime*(NumberOfBins-1)),
        hist2d_occupancy("occupancy", 99, 1.0, 100.0, 8, 1.0, 9.0)
{
	// Data
	ahdc = new AhdcDetector();
	adcCut = 0;
        cut_adcOffset_min = 0;
        cut_adcOffset_max = 4095;
        cut_leadingEdgeTime_min = 0;
        cut_leadingEdgeTime_max = 950;
        cut_timeOverThreshold_min = 0;
        cut_timeOverThreshold_max = 950;
        cut_timeMax_min = 0;
        cut_timeMax_max = 950;
	// Widgets
	set_title("ALERT monitoring");
	set_default_size(1600,990);
	set_child(VBox_main);
	
	/********************
	 * HEADER
	 * *****************/
	VBox_main.append(VBox_header);
	VBox_header.append(Label_header);
	Label_header.set_text("No file selected");
	Label_header.set_hexpand(true);
	/*******************
	 * BODY
	 * ****************/
	VBox_main.append(VBox_body);
	VBox_body.append(Book);
	Book.set_margin(10);
	Book.set_expand(true);
	Book.signal_switch_page().connect(sigc::mem_fun(*this,&Window::on_book_switch_page) );
	
	// Page 0
	Book.append_page(HBox_eventViewer, "Event Viewer");
	HBox_eventViewer.append(Grid_eventViewer);
	Grid_eventViewer.set_column_homogeneous(true);
	Grid_eventViewer.set_row_homogeneous(true);
	Grid_eventViewer.attach(DrawingArea_event,1,1);
	DrawingArea_event.set_draw_func(sigc::mem_fun(*this, &Window::on_draw_event) );
	Grid_eventViewer.attach(Grid_waveforms,2,1);
	Grid_waveforms.set_expand(true);
	Grid_waveforms.set_column_homogeneous(true);
	Grid_waveforms.set_row_homogeneous(true);
	// Page 1
	Book.append_page(Grid_waveformsPerLayer, "Waveforms Per Layer");
	Grid_waveformsPerLayer.set_expand(true);
	Grid_waveformsPerLayer.set_column_homogeneous(true);
	Grid_waveformsPerLayer.set_row_homogeneous(true);
	// Page 2
	Book.append_page(HBox_histograms, "Histograms");
	HBox_histograms.append(Grid_histograms);
	Grid_histograms.set_expand(true);
        Grid_histograms.set_column_homogeneous(true);
        Grid_histograms.set_row_homogeneous(true);
	hist1d_adcMax.set_xtitle("adc");
        hist1d_leadingEdgeTime.set_xtitle("time (ns)");
        hist1d_timeOverThreshold.set_xtitle("time (ns)");
        hist1d_timeMax.set_xtitle("time (ns)");
	hist1d_adcOffset.set_xtitle("adc");
        hist1d_constantFractionTime.set_xtitle("time (ns)");
	// page 3
	Book.append_page(Grid_occupancy, "Occupancy");
	Grid_occupancy.set_expand(true);
	Grid_occupancy.set_column_homogeneous(true);
	Grid_occupancy.set_row_homogeneous(true);
	TextBuffer_occupancy = Gtk::TextBuffer::create();
	// page 4
	Book.append_page(HBox_3Dview, "3D view");
	HBox_3Dview.append(VBox_3Dview_settings);
	//VBox_3Dview_settings.set_size_request(300,-1);
		// alpha
	auto HBox_Scale_alpha = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	HBox_Scale_alpha->set_margin_start(10);
	auto Label_alpha = Gtk::make_managed<Gtk::Label>();
	Label_alpha->set_markup("<b> alpha </b>");
	VBox_3Dview_settings.append(*HBox_Scale_alpha);
	HBox_Scale_alpha->append(*Label_alpha);
	HBox_Scale_alpha->append(Scale_alpha);
	Scale_alpha.set_hexpand(true);
	Scale_alpha.set_draw_value();
	Scale_alpha.set_digits(2);
	Adjustment_alpha->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_alpha->get_value();
				this->angle_alpha = val;
				DrawingArea_3Dview.queue_draw();
			});
		// beta
	auto HBox_Scale_beta = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	HBox_Scale_beta->set_margin_start(10);
	auto Label_beta = Gtk::make_managed<Gtk::Label>();
	Label_beta->set_markup("<b> beta </b>");
	VBox_3Dview_settings.append(*HBox_Scale_beta);
	HBox_Scale_beta->append(*Label_beta);
	HBox_Scale_beta->append(Scale_beta);
	Scale_beta.set_hexpand(true);
	Scale_beta.set_draw_value();
	Scale_beta.set_digits(2);
	Adjustment_beta->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_beta->get_value();
				this->angle_beta = val;
				DrawingArea_3Dview.queue_draw();
			});
		// gamma
	auto HBox_Scale_gamma = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	HBox_Scale_gamma->set_margin_start(10);
	auto Label_gamma = Gtk::make_managed<Gtk::Label>();
	Label_gamma->set_markup("<b> gamma </b>");
	VBox_3Dview_settings.append(*HBox_Scale_gamma);
	HBox_Scale_gamma->append(*Label_gamma);
	HBox_Scale_gamma->append(Scale_gamma);
	Scale_gamma.set_hexpand(true);
	Scale_gamma.set_draw_value();
	Scale_gamma.set_digits(2);
	Adjustment_gamma->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_gamma->get_value();
				this->angle_gamma = val;
				DrawingArea_3Dview.queue_draw();
			});
		// DrawingArea
	auto Separator_3Dview = Gtk::make_managed<Gtk::Separator>(Gtk::Orientation::VERTICAL);
	HBox_3Dview.append(*Separator_3Dview);
	HBox_3Dview.append(DrawingArea_3Dview);
	DrawingArea_3Dview.set_size_request(1400,-1);
	DrawingArea_3Dview.set_draw_func(sigc::mem_fun(*this, &Window::on_draw_3Dview) );
	// Page 4 (test)
	Book.append_page(DrawingArea_test, "Test");
	DrawingArea_test.set_draw_func(sigc::mem_fun(*this, &Window::on_draw_test) );

	/******************
	 * FOOTER
	 * ***************/
	VBox_main.append(VBox_footer);
	VBox_footer.append(HBox_footer);
	HBox_footer.set_margin(10);
	// settings
	HBox_footer.append(Button_settings);
	//Button_settings.set_hexpand(true);
	auto HBox_settings = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,15);
	Button_settings.set_child(*HBox_settings);
	img_settings.set("../img/icon_settings2.png");
	HBox_settings->append(img_settings );
	HBox_settings->append(*Gtk::make_managed<Gtk::Label>("settings") );
	Button_settings.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_settings_clicked) );
	// prev
	HBox_footer.append(Button_prev);
	//Button_prev.set_hexpand(true);
	auto HBox_prev = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,15);
	Button_prev.set_child(*HBox_prev);
	img_prev.set("../img/icon_prev_off.png");
	HBox_prev->append(img_prev );
	HBox_prev->append(*Gtk::make_managed<Gtk::Label>("prev") );
	Button_prev.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_prev_clicked) );

	// next
	HBox_footer.append(Button_next);
	auto HBox_next = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,15);
	Button_next.set_child(*HBox_next);
	img_next.set("../img/icon_next_off.png");
	HBox_next->append(img_next);
	HBox_next->append(*Gtk::make_managed<Gtk::Label>("next") );
	Button_next.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_next_clicked) );

	// pause
	HBox_footer.append(Button_pause);
	auto HBox_pause = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,20);
	Button_pause.set_child(*HBox_pause);
	img_pause.set("../img/icon_pause_off.png");
	HBox_pause->append(img_pause );
	HBox_pause->append(*Gtk::make_managed<Gtk::Label>("pause") );
	Button_pause.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_pause_clicked) );

	// run
	HBox_footer.append(Button_run);
	auto HBox_run = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,20);
	Button_run.set_child(*HBox_run);
	img_run.set("../img/icon_run_off.png");
	HBox_run->append(img_run);
	HBox_run->append(*Gtk::make_managed<Gtk::Label>("run") );
	Button_run.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_run_clicked) );

	// zpos
	HBox_footer.append(Button_zpos);
	auto HBox_zpos = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,20);
	Button_zpos.set_child(*HBox_zpos);
	img_zpos.set("../img/icon_zpos.png");
	HBox_zpos->append(img_zpos);
	HBox_zpos->append(*Gtk::make_managed<Gtk::Label>("zpos") );
	Button_zpos.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_zpos_clicked) );
	
	// middle
	HBox_footer.append(HBox_info);
	HBox_info.set_hexpand(true);
	HBox_info.append(Label_info);
	Label_info.set_text("No data");

	// hipo4
	HBox_footer.append(Button_hipo4);
	auto HBox_hipo4 = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,20);
	Button_hipo4.set_child(*HBox_hipo4);
	img_hipo4.set("../img/icon_file_on.png");
	HBox_hipo4->append(img_hipo4);
	HBox_hipo4->append(*Gtk::make_managed<Gtk::Label>("hipo4") );
	Button_hipo4.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_hipo4_clicked) );

	// reset
	HBox_footer.append(Button_reset);
	auto HBox_reset = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,20);
	Button_reset.set_child(*HBox_reset);
	img_reset.set("../img/icon_reset_off.png");
	HBox_reset->append(img_reset);
	HBox_reset->append(*Gtk::make_managed<Gtk::Label>("reset") );
	Button_reset.signal_clicked().connect( sigc::mem_fun(*this, &Window::on_button_reset_clicked) );
	
	//VBox_footer.append(*Gtk::make_managed<Gtk::Label>("Footer") );
	
	/**************************
	 * Go back in eventViewer
	 * **********************/
	auto mouse_click = Gtk::GestureClick::create();
	mouse_click->signal_pressed().connect(sigc::mem_fun(*this, &Window::on_mouse_clicked));
	DrawingArea_event.add_controller(mouse_click);
}

/** Destructor */
Window::~Window() {
	// nothing
	delete ahdc;
}

void Window::on_button_settings_clicked(){
	std::cout << "Open settings tab ..." << std::endl;
	// Initialization
	auto Window_settings = Gtk::make_managed<Gtk::Window>();
	auto VBox_settings = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::VERTICAL);
	auto HBox_Scale_adcMax = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_adcOffset_min = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_adcOffset_max = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_leadingEdgeTime_min = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_leadingEdgeTime_max = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_timeOverThreshold_min = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_timeOverThreshold_max = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_timeMax_min = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_Scale_cut_timeMax_max = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	auto HBox_CheckButtons = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	// Window for settings
	Window_settings->set_title("Settings");
	Window_settings->set_default_size(500,300);
	Window_settings->set_child(*VBox_settings);
	// Scale_adcMax
	auto Label0 = Gtk::make_managed<Gtk::Label>();
	Label0->set_markup("<b> adcMax </b>");
	VBox_settings->append(*Label0);
	VBox_settings->append(*HBox_Scale_adcMax);
	HBox_Scale_adcMax->set_margin_start(10);
	HBox_Scale_adcMax->append(*Gtk::make_managed<Gtk::Label>("MIN"));
	HBox_Scale_adcMax->append(Scale_adcMax);
	Scale_adcMax.set_hexpand(true);
	Scale_adcMax.set_draw_value();
	Scale_adcMax.set_digits(0);
	Adjustment_adcMax->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_adcMax->get_value();
				this->adcCut = val;
			});
	auto Separator0 = Gtk::make_managed<Gtk::Separator>();
	VBox_settings->append(*Separator0);	
	// Scale_cut_adcOffset_min
	auto Label1 = Gtk::make_managed<Gtk::Label>();
	Label1->set_markup("<b> adcOffset </b>");
	VBox_settings->append(*Label1);
	VBox_settings->append(*HBox_Scale_cut_adcOffset_min);
        HBox_Scale_cut_adcOffset_min->set_margin_start(10);
        HBox_Scale_cut_adcOffset_min->append(*Gtk::make_managed<Gtk::Label>("MIN"));
        HBox_Scale_cut_adcOffset_min->append(Scale_cut_adcOffset_min);	
	Scale_cut_adcOffset_min.set_hexpand(true);
	Scale_cut_adcOffset_min.set_draw_value();
	Scale_cut_adcOffset_min.set_digits(0);
	Adjustment_cut_adcOffset_min->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_adcOffset_min->get_value();
				this->cut_adcOffset_min = val;
			});
	// Scale_cut_adcOffset_max
	VBox_settings->append(*HBox_Scale_cut_adcOffset_max);
        HBox_Scale_cut_adcOffset_max->set_margin_start(10);
        HBox_Scale_cut_adcOffset_max->append(*Gtk::make_managed<Gtk::Label>("MAX"));
        HBox_Scale_cut_adcOffset_max->append(Scale_cut_adcOffset_max);
	Scale_cut_adcOffset_max.set_hexpand(true);
	Scale_cut_adcOffset_max.set_draw_value();
	Scale_cut_adcOffset_max.set_digits(0);
	Adjustment_cut_adcOffset_max->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_adcOffset_max->get_value();
				this->cut_adcOffset_max = val;
			});
	auto Separator1 = Gtk::make_managed<Gtk::Separator>();
	VBox_settings->append(*Separator1);	
	// Scale_cut_leadingEdgeTime_min
	auto Label2 = Gtk::make_managed<Gtk::Label>();
	Label2->set_markup("<b> leadingEdgeTime </b>");
	VBox_settings->append(*Label2);
	VBox_settings->append(*HBox_Scale_cut_leadingEdgeTime_min);
        HBox_Scale_cut_leadingEdgeTime_min->set_margin_start(10);
        HBox_Scale_cut_leadingEdgeTime_min->append(*Gtk::make_managed<Gtk::Label>("MIN"));
        HBox_Scale_cut_leadingEdgeTime_min->append(Scale_cut_leadingEdgeTime_min);	
	Scale_cut_leadingEdgeTime_min.set_hexpand(true);
	Scale_cut_leadingEdgeTime_min.set_draw_value();
	Scale_cut_leadingEdgeTime_min.set_digits(0);
	Adjustment_cut_leadingEdgeTime_min->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_leadingEdgeTime_min->get_value();
				this->cut_leadingEdgeTime_min = val;
			});
	// Scale_cut_leadingEdgeTime_max
	VBox_settings->append(*HBox_Scale_cut_leadingEdgeTime_max);
        HBox_Scale_cut_leadingEdgeTime_max->set_margin_start(10);
        HBox_Scale_cut_leadingEdgeTime_max->append(*Gtk::make_managed<Gtk::Label>("MAX"));
        HBox_Scale_cut_leadingEdgeTime_max->append(Scale_cut_leadingEdgeTime_max);
	Scale_cut_leadingEdgeTime_max.set_hexpand(true);
	Scale_cut_leadingEdgeTime_max.set_draw_value();
	Scale_cut_leadingEdgeTime_max.set_digits(0);
	Adjustment_cut_leadingEdgeTime_max->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_leadingEdgeTime_max->get_value();
				this->cut_leadingEdgeTime_max = val;
			});
	auto Separator2 = Gtk::make_managed<Gtk::Separator>();
	VBox_settings->append(*Separator2);	
	// Scale_cut_timeOverThreshold_min
	auto Label3 = Gtk::make_managed<Gtk::Label>();
	Label3->set_markup("<b> timeOverThreshold </b>");
        VBox_settings->append(*Label3);
	VBox_settings->append(*HBox_Scale_cut_timeOverThreshold_min);
        HBox_Scale_cut_timeOverThreshold_min->set_margin_start(10);
        HBox_Scale_cut_timeOverThreshold_min->append(*Gtk::make_managed<Gtk::Label>("MIN"));
        HBox_Scale_cut_timeOverThreshold_min->append(Scale_cut_timeOverThreshold_min);
	Scale_cut_timeOverThreshold_min.set_hexpand(true);
	Scale_cut_timeOverThreshold_min.set_draw_value();
	Scale_cut_timeOverThreshold_min.set_digits(0);
	Adjustment_cut_timeOverThreshold_min->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_timeOverThreshold_min->get_value();
				this->cut_timeOverThreshold_min = val;
			});
	// Scale_cut_timeOverThreshold_max
	VBox_settings->append(*HBox_Scale_cut_timeOverThreshold_max);
        HBox_Scale_cut_timeOverThreshold_max->set_margin_start(10);
        HBox_Scale_cut_timeOverThreshold_max->append(*Gtk::make_managed<Gtk::Label>("MAX"));
        HBox_Scale_cut_timeOverThreshold_max->append(Scale_cut_timeOverThreshold_max);
	Scale_cut_timeOverThreshold_max.set_hexpand(true);
	Scale_cut_timeOverThreshold_max.set_draw_value();
	Scale_cut_timeOverThreshold_max.set_digits(0);
	Adjustment_cut_timeOverThreshold_max->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_timeOverThreshold_max->get_value();
				this->cut_timeOverThreshold_max = val;
			});
	auto Separator3 = Gtk::make_managed<Gtk::Separator>();
	VBox_settings->append(*Separator3);	
	// Scale_cut_timeMax_min
	auto Label4 = Gtk::make_managed<Gtk::Label>();
	Label4->set_markup("<b> timeMax </b>");
        VBox_settings->append(*Label4);
	VBox_settings->append(*HBox_Scale_cut_timeMax_min);
        HBox_Scale_cut_timeMax_min->set_margin_start(10);
        HBox_Scale_cut_timeMax_min->append(*Gtk::make_managed<Gtk::Label>("MIN"));
        HBox_Scale_cut_timeMax_min->append(Scale_cut_timeMax_min);
	Scale_cut_timeMax_min.set_hexpand(true);
	Scale_cut_timeMax_min.set_draw_value();
	Scale_cut_timeMax_min.set_digits(0);
	Adjustment_cut_timeMax_min->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_timeMax_min->get_value();
				this->cut_timeMax_min = val;
			});
	// Scale_cut_timeMax_max
	VBox_settings->append(*HBox_Scale_cut_timeMax_max);
        HBox_Scale_cut_timeMax_max->set_margin_start(10);
        HBox_Scale_cut_timeMax_max->append(*Gtk::make_managed<Gtk::Label>("MAX"));
        HBox_Scale_cut_timeMax_max->append(Scale_cut_timeMax_max);
	Scale_cut_timeMax_max.set_hexpand(true);
	Scale_cut_timeMax_max.set_draw_value();
	Scale_cut_timeMax_max.set_digits(0);
	Adjustment_cut_timeMax_max->signal_value_changed().connect([this] () -> void {
				const double val = this->Adjustment_cut_timeMax_max->get_value();
				this->cut_timeMax_max = val;
			});
	// A lot of CheckButton(s)
	auto Separator4 = Gtk::make_managed<Gtk::Separator>();
	VBox_settings->append(*Separator4);
	VBox_settings->append(*HBox_CheckButtons);
	HBox_CheckButtons->set_margin_start(10);
	HBox_CheckButtons->set_margin_bottom(10);
	HBox_CheckButtons->append(*Gtk::make_managed<Gtk::Label>("Mask wires"));
	HBox_CheckButtons->append(CheckButton_mask_wires);
	CheckButton_mask_wires.signal_toggled().connect([this] () -> void {
				flag_mask_wires = this->CheckButton_mask_wires.get_active();
				if (flag_mask_wires) {
					printf("Flag_mask_wires is activated\n");
					DrawingArea_event.queue_draw();
				} else {
					printf("Flag_mask_wires is desactivated\n");
					DrawingArea_event.queue_draw();
				}
			});
	Window_settings->show();
}

void Window::on_button_prev_clicked(){
	std::cout << "Process prev event ..." << std::endl;
	// update icons
	//img_next.set("../img/icon_prev_off.png"); img_prev.queue_draw();
	//img_next.set("../img/icon_next_off.png"); img_next.queue_draw();
	//img_pause.set("../img/icon_pause_off.png"); img_pause.queue_draw();
	//img_run.set("../img/icon_run_off.png"); img_run.queue_draw();
	//img_hipo4.set("../img/icon_file_off.png"); img_hipo4.queue_draw();
	//img_reset.set("../img/icon_reset_off.png"); img_reset.queue_draw();
}

void Window::on_button_next_clicked(){
	std::cout << "Process next event ..." << std::endl;
	is_paused = false;
	is_reset = false;
	// update icons
	img_next.set("../img/icon_prev_off.png"); img_prev.queue_draw();
	img_next.set("../img/icon_next_on.png"); img_next.queue_draw();
	img_pause.set("../img/icon_pause_on.png"); img_pause.queue_draw();
	img_run.set("../img/icon_run_on.png"); img_run.queue_draw();
	img_hipo4.set("../img/icon_file_off.png"); img_hipo4.queue_draw();
	img_reset.set("../img/icon_reset_on.png"); img_reset.queue_draw();
	if (filename.size() == 0) {
		return;
	}
	if (hipo_nEvent == 0) {
		hipo_reader.open(filename.c_str());
		hipo_reader.readDictionary(hipo_factory);
		adcBank = hipo::bank(hipo_factory.getSchema("AHDC::adc"));
		wfBank = hipo::bank(hipo_factory.getSchema("AHDC::wf"));
		trackBank = hipo::bank(hipo_factory.getSchema("AHDC::kftrack"));
		hipo_nEventMax = hipo_reader.getEntries();
	}
	Glib::signal_timeout().connect([this] () -> bool {
				if (is_paused) {return false;}
				if (this->dataEventAction()) {
					if (nWF == 0) {return true;} // continue the timeout
					else {return false;} // stop the timeout
				}
				else {return false;} // stop the timeout
			}, 20); // call every 10 ms
}

void Window::on_button_pause_clicked(){
	std::cout << "Pause ..." << std::endl;
	is_paused = true;
	// update icons
	//img_prev.set("../img/icon_prev_off.png"); img_prev.queue_draw();
	//img_next.set("../img/icon_next_off.png"); img_next.queue_draw();
	img_pause.set("../img/icon_pause_off.png"); img_pause.queue_draw();
	//img_run.set("../img/icon_run_off.png"); img_run.queue_draw();
	//img_hipo4.set("../img/icon_file_off.png"); img_hipo4.queue_draw();
	//img_reset.set("../img/icon_reset_off.png"); img_reset.queue_draw();
	drawWaveforms();
	drawWaveformsPerLayer();
}

void Window::on_button_run_clicked(){
	std::cout << "Run ..." << std::endl;
	is_paused = false;
	is_reset = false;
	// update icons
	img_prev.set("../img/icon_prev_off.png"); img_prev.queue_draw();
	img_next.set("../img/icon_next_on.png"); img_next.queue_draw();
	img_pause.set("../img/icon_pause_on.png"); img_pause.queue_draw();
	img_run.set("../img/icon_run_on.png"); img_run.queue_draw();
	img_hipo4.set("../img/icon_file_off.png"); img_hipo4.queue_draw();
	img_reset.set("../img/icon_reset_on.png"); img_reset.queue_draw();
	if (filename.size() == 0) {
                return;
        }
        if (hipo_nEvent == 0) {
                hipo_reader.open(filename.c_str());
		hipo_reader.readDictionary(hipo_factory);
		adcBank = hipo::bank(hipo_factory.getSchema("AHDC::adc"));
		wfBank = hipo::bank(hipo_factory.getSchema("AHDC::wf"));
		trackBank = hipo::bank(hipo_factory.getSchema("AHDC::kftrack"));
		hipo_nEventMax = hipo_reader.getEntries();
        }
	Glib::signal_timeout().connect([this] () -> bool {
				if (is_paused || is_reset) {return false;}
				if (this->dataEventAction()) {
					return true; // continue the timeout
					//return true && (hipo_nEvent < 10000); 
				}
				else {return false;} // stop the timeout
                        }, 20); // call every 10 ms
}

void Window::on_button_zpos_clicked(){
	auto Window_zpos = Gtk::make_managed<Gtk::Window>();
	auto HBox_zpos = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL,10);
	HBox_zpos->set_margin(10);
	Window_zpos->set_title("zpos");
	Window_zpos->set_default_size(500,30);
	Window_zpos->set_child(*HBox_zpos);
	HBox_zpos->append(*Gtk::make_managed<Gtk::Label>("zpos"));
	HBox_zpos->append(Scale_zpos);
	Scale_zpos.set_hexpand(true);
	Scale_zpos.set_draw_value();
	Scale_zpos.set_digits(0);
	Adjustment_zpos->signal_value_changed().connect(sigc::mem_fun(*this, &Window::on_zpos_value_changed));
	Window_zpos->show();
}

void Window::on_button_hipo4_clicked(){
	std::cout << "Open file dialog ..." << std::endl;

	auto dialog = Gtk::FileDialog::create();

	// Add filters, so that only certain file types can be selected:
	auto filters = Gio::ListStore<Gtk::FileFilter>::create();

	auto filter_hipo4 = Gtk::FileFilter::create();
	filter_hipo4->set_name("Hipo files");
	filter_hipo4->add_pattern("*.hipo");
	filters->append(filter_hipo4);
	
	auto filter_evio = Gtk::FileFilter::create();
	filter_evio->set_name("Evio files");
	filter_evio->add_pattern("*.evio");
	filters->append(filter_evio);
	
	auto filter_text = Gtk::FileFilter::create();
	filter_text->set_name("Text files");
	filter_text->add_mime_type("text/plain");
	filters->append(filter_text);

	auto filter_cpp = Gtk::FileFilter::create();
	filter_cpp->set_name("C/C++ files");
	filter_cpp->add_mime_type("text/x-c");
	filter_cpp->add_mime_type("text/x-c++");
	filter_cpp->add_mime_type("text/x-c-header");
	filters->append(filter_cpp);

	auto filter_any = Gtk::FileFilter::create();
	filter_any->set_name("Any files");
	filter_any->add_pattern("*");
	filters->append(filter_any);

	dialog->set_filters(filters);

	// Show the dialog and wait for a user response:
	dialog->open(sigc::bind(sigc::mem_fun(*this, &Window::on_file_dialog_finish), dialog));
}
void Window::on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog) {
	// Handle the response:
	try
	{
		auto file = dialog->open_finish(result);

		// Notice that this is a std::string, not a Glib::ustring.
		filename = file->get_path();
		std::cout << "File selected : " <<  filename << std::endl;
		
		// Possible actions
		img_prev.set("../img/icon_prev_off.png"); img_prev.queue_draw();
		img_next.set("../img/icon_next_on.png"); img_next.queue_draw();
		img_pause.set("../img/icon_pause_off.png"); img_pause.queue_draw();
		img_run.set("../img/icon_run_on.png"); img_run.queue_draw();
		img_hipo4.set("../img/icon_file_off.png"); img_hipo4.queue_draw();
		img_reset.set("../img/icon_reset_off.png"); img_reset.queue_draw();
		// label info
		this->Label_header.set_text(TString::Format("File selected : %s", filename.c_str()).Data() );
		Label_info.queue_draw();
	}
	catch (const Gtk::DialogError& err)
	{
		// Can be thrown by dialog->open_finish(result).
		std::cout << "No file selected. " << err.what() << std::endl;
	}
	catch (const Glib::Error& err)
	{
		std::cout << "Unexpected exception. " << err.what() << std::endl;
	}
}

void Window::on_button_reset_clicked(){
	std::cout << "Reset ..." << std::endl;
	is_reset = true;
	filename = "";
	// update icons
	img_prev.set("../img/icon_prev_off.png"); img_prev.queue_draw();
	img_next.set("../img/icon_next_off.png"); img_next.queue_draw();
	img_pause.set("../img/icon_pause_off.png"); img_pause.queue_draw();
	img_run.set("../img/icon_run_off.png"); img_run.queue_draw();
	img_hipo4.set("../img/icon_file_on.png"); img_hipo4.queue_draw();
	img_reset.set("../img/icon_reset_off.png"); img_reset.queue_draw();
	// ...
	hipo_nEvent = 0;
	hipo_nEventMax = 1;
	clearAhdcData();
	ListOfAdc.clear();
	// Reset eventViewer waveforms
	DrawingArea_event.queue_draw();
	drawWaveforms();
	drawWaveformsPerLayer();
	// Reset histograms
	hist1d_adcMax.reset();
	hist1d_leadingEdgeTime.reset();
	hist1d_timeOverThreshold.reset();
	hist1d_timeMax.reset();
	hist1d_adcOffset.reset();
        hist1d_constantFractionTime.reset();
	hist2d_occupancy.reset();
	drawHistograms();
	// Reset occupancy
	drawOccupancy();
	// Update Label
	Label_info.set_text("No data"); Label_info.queue_draw();
	Label_header.set_text("No file selected"); Label_header.queue_draw();
}

void Window::on_book_switch_page(Gtk::Widget * page, guint page_num) { 
	std::string page_name;
	switch (page_num) {
		case 0 :
			page_name = "Event Viewer";
			DrawingArea_event.queue_draw();
			drawWaveforms();
			break;
		case 1 :
			page_name = "Waveforms Per Layer";
			drawWaveformsPerLayer();
			break;
		case 2 :
			page_name = "Histograms";
			drawHistograms();
			break;
		case 3 :
			page_name = "Occupancy";
			drawOccupancy();
			break;
		case 4 :
			page_name = "3D view";
			break;
		default :
			page_name = "Unknown process";
	}
	std::cout << "Switch to " << page_name << " tab ..." << std::endl;
}

void Window::on_draw_event(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height){
	fCanvas canvas(width, height, -80, 80, -80, 80);
	int window_size = std::min(width,height);
	canvas.set_top_margin(0.05*window_size);
	canvas.set_bottom_margin(0.05*window_size);
	canvas.set_left_margin(0.05*window_size);
	//canvas.set_right_margin(0.05*window_size);
	canvas.set_right_margin(0.15*window_size);
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
	char face_name[50];
	sprintf(face_name, "(x,y) view in plane z = %.2lf mm", zpos);
	canvas.draw_title(cr, face_name);
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
	// Draw axis orientation
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->set_line_width(0.003*seff);
		// x axis
	cr->move_to(x2w(-70), y2h(-70));
	cr->line_to(x2w(-70+20), y2h(-70));
	cr->stroke();
	cr->move_to(x2w(-70+20), y2h(-70));
	cr->line_to(x2w(-70+20) + 20*cos(M_PI-0.25*M_PI), y2h(-70) + 10*sin(M_PI-0.25*M_PI));
	cr->line_to(x2w(-70+20) + 20*cos(M_PI+0.25*M_PI), y2h(-70) + 10*sin(M_PI+0.25*M_PI));
	cr->close_path();
	cr->fill();
	{ // x axis label
		cr->set_source_rgb(0.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.2*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("x", te);
		cr->move_to(x2w(-70+20)+3, y2h(-70)+3);
		cr->show_text("x");
	}
		// y axis
	cr->move_to(x2w(-70), y2h(-70));
	cr->line_to(x2w(-70), y2h(-70+20));
	cr->stroke();
	cr->move_to(x2w(-70), y2h(-70+20));
	cr->line_to(x2w(-70) + 10*cos(M_PI+0.25*M_PI), y2h(-70+20) - 20*sin(M_PI+0.25*M_PI)); // minus because of the orientation :-(
	cr->line_to(x2w(-70) + 10*cos(-0.25*M_PI), y2h(-70+20) - 20*sin(-0.25*M_PI));
	cr->close_path();
	cr->fill();
	{ // y axis label
		cr->set_source_rgb(0.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.2*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("y", te);
		cr->move_to(x2w(-70)-3-te.width, y2h(-70+20)+3);
		cr->show_text("y");
	}
		// z axis
	cr->set_source_rgb(1.0, 0.0, 0.0);
	cr->move_to(x2w(-70)+10, y2h(-70));
	cr->arc(x2w(-70),y2h(-70), 10, 0, 2*M_PI);
	cr->stroke();
	cr->move_to(x2w(-70)+3, y2h(-70));
	cr->arc(x2w(-70),y2h(-70), 3, 0, 2*M_PI);
	cr->fill();
	{ // z axis label
		cr->set_source_rgb(1.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.5*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("z", te);
		cr->move_to(x2w(-70) + 13, y2h(-70) + 13 + te.height);
		cr->show_text("z");
	}
		
	// Draw color palette
	double zmin = 0.0;
	double zmax = 4095.0;
	if (ListOfAdc.size() > 0) {
		zmin = ListOfAdc[0];
		zmax = zmin;
	}
	for (double adc : ListOfAdc) {
		zmin = (zmin < adc) ? zmin : adc;
		zmax = (zmax > adc) ? zmax : adc;	
	}
	//printf("zmin : %.0lf, zmax : %.0lf\n", zmin, zmax);
	fColorPalette Palette(5, 1);
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->set_line_width(0.005*seff);
	cr->move_to(x2w(80) + 0.01*window_size, y2h(-80)); // 0.1*window_size is the extra margin of right_margin
	cr->line_to(x2w(80) + 0.1*window_size, y2h(-80));
	cr->line_to(x2w(80) + 0.1*window_size, y2h(80));
	cr->stroke();
	for (int i = 0; i < Palette.get_ncolors(); i++) {
		int ncolors = Palette.get_ncolors();
		fColor rgb = Palette.get_color(i);
		double slope = 160.0/(ncolors); // == 80 - (-80)
		double z = slope*(i-0.0) - 80.0;
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->move_to(x2w(80) + 0.01*window_size, y2h(z));
		cr->line_to(x2w(80) +  0.1*window_size, y2h(z));
		cr->line_to(x2w(80) +  0.1*window_size, y2h(z+1.0*slope));
		cr->line_to(x2w(80) + 0.01*window_size, y2h(z+1.0*slope));
		cr->close_path();
		cr->fill();
	}
	{ // Draw zmax
		char buffer[50];
		sprintf(buffer, "%.0lf", zmax);
		//fColor rgb = Palette.get_color(Palette.get_ncolors()-1);
		fColor rgb({0.0,0.0,0.0});
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents(buffer, te);
		cr->move_to(x2w(80) + 0.05*window_size - 0.5*te.width, y2h(80) - 0.1*canvas.get_top_margin());
		cr->show_text(buffer);
	}
	// Draw zmin
	{ 
		char buffer[50];
		sprintf(buffer, "%.0lf", zmin);
		fColor rgb = Palette.get_color(0);
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents(buffer, te);
		cr->move_to(x2w(80) + 0.05*window_size - 0.5*te.width, y2h(-80) + 0.1*canvas.get_top_margin() + te.height);
		cr->show_text(buffer);
	}	
	// Draw AHDC geometry	
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
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
					//printf("%lf %lf %lf\n", wire->x, wire->y, wire->z);
					cr->move_to(x2w(wire->x) + marker_size, y2h(wire->y)); // correspond to angle == 0 
					cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
					cr->stroke();
					//if ((wire->x == 0) && (wire->y < 0)) {
					if (w == 0) {
						// these wires have a component id == 1 (it is the start of the numerotation)
						cr->set_line_width(0.002*seff);
						cr->set_source_rgba(0.0, 0.0, 0.2, 1.0);
						cr->move_to(x2w(wire->x) + marker_size, y2h(wire->y));
						cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
						cr->stroke();
					}
					// show activated wires
					double adcMax            = wire->pulse.get_adcMax();
					if (adcMax > 0) {
						bool flag = true;
						double adcOffset         = wire->pulse.get_adcOffset();
						double leadingEdgeTime   = wire->pulse.get_leadingEdgeTime();
						double timeOverThreshold = wire->pulse.get_timeOverThreshold();
						double timeMax           = wire->pulse.get_timeMax();
						if (flag_mask_wires) {
							flag = (adcMax >= adcCut) &&
								(adcOffset >= cut_adcOffset_min) && 
								(adcOffset <= cut_adcOffset_max) && 
								(leadingEdgeTime >= cut_leadingEdgeTime_min) &&
								(leadingEdgeTime <= cut_leadingEdgeTime_max) &&
								(timeOverThreshold >= cut_timeOverThreshold_min) && 
								(timeOverThreshold <= cut_timeOverThreshold_max) &&
								(timeMax >= cut_timeMax_min) &&
								(timeMax <= cut_timeMax_max); 
						}
						if (!flag) continue;
						if (zmin != zmax) {
							int ncolors = Palette.get_ncolors();
							double slope = (1.0*(ncolors-1))/(zmax-zmin);
							int color = floor(slope*(adcMax-zmin) + 0.0);
							//printf("zmin : %.0lf, zmax : %.0lf\n", zmin, zmax);
							//printf("color : %d\n", color);
							fColor rgb = Palette.get_color(color);
							cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
						} else {
							fColor rgb = Palette.get_color(0);
							cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
						}
						double marker_size = std::min(2.0*weff/(x_end-x_start), 2.0*heff/(y_end-y_start));
						cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
						cr->fill();
					}
				}
			}
		}
	}
}

void Window::on_draw_test(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height){
	cr->save();
	fColorPalette Palette(5, 2);	
	int ncolors = Palette.get_ncolors();
	fCanvas canvas(width, height, 0.0, 1.0*ncolors, 0.0, 1.0);
	canvas.define_coord_system(cr);
	// x coord to width
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	// y coord to height
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};
	for (int i = 0; i < ncolors; i++) {
		fColor rgb = Palette.get_color(i);
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->set_line_width(0.01*canvas.get_seff());
		cr->move_to(x2w(1.0*i), y2h(0.0));
		cr->line_to(x2w(1.0*(i+1)), y2h(0.0));
		cr->line_to(x2w(1.0*(i+1)), y2h(1.0));
		cr->line_to(x2w(1.0*i), y2h(1.0));
		cr->close_path();
		//cr->stroke_preserve();
		cr->fill();
	}	
	canvas.do_not_draw_secondary_stick();
	canvas.draw_frame(cr);
	//canvas.set_title_size(0.6);
	canvas.draw_title(cr, "Example (test fColorPalette)");
	canvas.draw_xtitle(cr, "color");
	canvas.draw_ytitle(cr, "");
	cr->restore();
}

void Window::cairo_plot_waveform(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height, AhdcWire* wire, std::string annotation){
	// Determine the min and the max of the data	
	std::vector<double> vy = wire->pulse.get_samples();
	if (NumberOfBins != (int) vy.size()) {
		return ;
	}
	std::vector<double> vx;
	for (int i = 0; i < NumberOfBins; i++){
		vx.push_back(samplingTime*i);
	}
	double xmin = vx[0], xmax = vx[0];
	double ymin = vy[0], ymax = vy[0];
	for (int i = 0; i < NumberOfBins; i++){
		xmin = (xmin < vx[i]) ? xmin : vx[i];
 		xmax = (xmax > vx[i]) ? xmax : vx[i];
		ymin = (ymin < vy[i]) ? ymin : vy[i];
		ymax = (ymax > vy[i]) ? ymax : vy[i];
	}
	fCanvas canvas(width, height, xmin, xmax, ymin, ymax);
	canvas.define_coord_system(cr);
	canvas.draw_title(cr, "");
	canvas.draw_xtitle(cr, "time (ns)");
	canvas.draw_ytitle(cr, "adc");
	// x coord to width
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	// y coord to height
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};

	// draw frame and axis
	canvas.draw_frame(cr); // can be placed here or at the end, the most important is to have defined `canvas.define_coord_system(cr);` before draw plots
	
	int seff = canvas.get_seff();
	int heff = canvas.get_heff();
	int weff = canvas.get_weff();	
	// Draw points
	cr->set_source_rgb(0.0, 0.0, 1.0);
	cr->set_line_width(0.01*seff);
	cr->move_to(x2w(vx[0]),y2h(vy[0]));
	for (int i = 1; i < NumberOfBins; i++) {
		// draw a line between points i and i-1
		cr->line_to(x2w(vx[i]),y2h(vy[i]));
	}
	cr->stroke();
	
	// ___________________________
	// Show decoded values
	// ___________________________
	//int binOffset = wire->pulse.get_binOffset();
	//double timeMax = wire->pulse.get_timeMax();
	double leadingEdgeTime = wire->pulse.get_leadingEdgeTime();
	double constantFractionTime = wire->pulse.get_constantFractionTime();
	double timeOverThreshold = wire->pulse.get_timeOverThreshold();
	double adcMax = wire->pulse.get_adcMax();
	double adcOffset = wire->pulse.get_adcOffset();
	
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
	cr->move_to(x2w(leadingEdgeTime), y2h(adcOffset + adcMax*amplitudeFractionCFA));
	cr->line_to(x2w(leadingEdgeTime + timeOverThreshold), y2h(adcOffset + adcMax*amplitudeFractionCFA));
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
	cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
	cr->set_font_size(seff*0.1);
	cr->move_to(weff*0.7, -heff*0.8);
	cr->show_text(annotation);
	
}

bool Window::dataEventAction() {
	// hipo4
	if (hipo_reader.next()) {
		hipo_reader.read(hipo_event);
		hipo_event.getStructure(adcBank);
		hipo_event.getStructure(wfBank);
		hipo_event.getStructure(trackBank);
		hipo_nEvent++;
		// loop over hits
		clearAhdcData();
		ListOfAdc.clear();
		bool doIshowWF = false;
		nWF = 0;
		for (int col = 0; col < wfBank.getRows(); col++){
			//int sector = wfBank.getInt("sector", col);	
			int layer = wfBank.getInt("layer", col);
			int component = wfBank.getInt("component", col);
			std::vector<double> samples;
			for (int bin=0; bin < NumberOfBins; bin++){
				std::string binName = "s" + std::__cxx11::to_string(bin+1);
				short value = wfBank.getInt(binName.c_str(), col);
				samples.push_back(value);
			}
			// fill histograms
			double factor = 1.0; // if == samplingTime, (the unit in bin number)
			double timeMax = this->adcBank.getFloat("time", col)/factor;
                        double leadingEdgeTime = this->adcBank.getFloat("leadingEdgeTime", col)/factor;
                        double timeOverThreshold = this->adcBank.getFloat("timeOverThreshold", col)/factor;
                        double constantFractionTime = this->adcBank.getFloat("constantFractionTime", col)/factor;
                        double adcMax = this->adcBank.getInt("ADC", col); // expected adcMax without adcOffset
                        double adcOffset = this->adcBank.getInt("ped", col);
                        //double integral = this->adcBank.getInt("integral", col);
			bool flag = (adcMax >= adcCut) &&
				(adcOffset >= cut_adcOffset_min) && 
				(adcOffset <= cut_adcOffset_max) && 
				(leadingEdgeTime >= cut_leadingEdgeTime_min) &&
				(leadingEdgeTime <= cut_leadingEdgeTime_max) &&
				(timeOverThreshold >= cut_timeOverThreshold_min) && 
				(timeOverThreshold <= cut_timeOverThreshold_max) &&
				(timeMax >= cut_timeMax_min) &&
				(timeMax <= cut_timeMax_max);
                       	//bool flag = (adcMax >= adcCut) && (adcOffset >= cut_adcOffset_min) && (adcOffset <= cut_adcOffset_max) && (leadingEdgeTime >= cut_leadingEdgeTime_min) && (leadingEdgeTime <= cut_leadingEdgeTime_max) && (timeOverThreshold >= cut_timeOverThreshold_min) && (timeOverThreshold <= cut_timeOverThreshold_max) && (timeMax >= cut_timeMax_min) && (timeMax <= cut_timeMax_max); 
			if (flag) {
					hist1d_adcMax.fill(adcMax);
					hist1d_leadingEdgeTime.fill(leadingEdgeTime);
					hist1d_timeOverThreshold.fill(timeOverThreshold);
					hist1d_timeMax.fill(timeMax);
					hist1d_adcOffset.fill(adcOffset);
					hist1d_constantFractionTime.fill(constantFractionTime);
					hist2d_occupancy.fill(component, layer2number(layer));
					doIshowWF = true; // at least one wire has reach the adc cut; so show the whole event (see next loop for)
					//doIshowWF = true && (layer != 51) && (layer !=42); // prevent the last layer to trigger
			}
		}
		// Criteria to show all the event 	
		if (doIshowWF) {
			for (int col = 0; col < wfBank.getRows(); col++){
				int sector = wfBank.getInt("sector", col);	
				int layer = wfBank.getInt("layer", col);
				int component = wfBank.getInt("component", col);
				int binOffset = wfBank.getInt("time", col);
				std::vector<double> samples;
				for (int bin=0; bin < binOffset; bin++){
					samples.push_back(0.0);
				}
				for (int bin=0; bin < NumberOfBins - binOffset; bin++){
					std::string binName = "s" + std::__cxx11::to_string(bin+1);
					short value = wfBank.getInt(binName.c_str(), col);
					samples.push_back(value);
				}
				double factor = 1.0; // if == samplingTime, (the unit in bin number)
				double timeMax = this->adcBank.getFloat("time", col)/factor;
				double leadingEdgeTime = this->adcBank.getFloat("leadingEdgeTime", col)/factor;
				double timeOverThreshold = this->adcBank.getFloat("timeOverThreshold", col)/factor;
				double constantFractionTime = this->adcBank.getFloat("constantFractionTime", col)/factor;
				double adcMax = this->adcBank.getInt("ADC", col); // expected adcMax without adcOffset
				double adcOffset = this->adcBank.getInt("ped", col);
				double integral = this->adcBank.getInt("integral", col);
				AhdcWire *wire = ahdc->GetSector(sector-1)->GetSuperLayer((layer/10)-1)->GetLayer((layer%10)-1)->GetWire(component-1);	
				wire->pulse.set_adcMax(adcMax);
				wire->pulse.set_adcOffset(adcOffset);
				wire->pulse.set_integral(integral);
				wire->pulse.set_timeMax(timeMax);
				wire->pulse.set_leadingEdgeTime(leadingEdgeTime);
				wire->pulse.set_timeOverThreshold(timeOverThreshold);
				wire->pulse.set_constantFractionTime(constantFractionTime);
				wire->pulse.set_binOffset(binOffset);
				wire->pulse.set_samples(samples);
				wire->occ += 1;
				ListOfAdc.push_back(adcMax);
				nWF++;
			}
		}
		// Update drawings
		int tab_number = Book.get_current_page();
		if      (tab_number == 0) {
			DrawingArea_event.queue_draw();
			drawWaveforms();
		}
		else if (tab_number == 1) {
			drawWaveformsPerLayer();
		}
		else if (tab_number == 2) {
			drawHistograms();
		}
		else if (tab_number == 3) {
			drawOccupancy();
		}
		else {
			// do nothing
		}
		Label_info.set_text(TString::Format("Progress : %.2lf %%, Event number : %lu/%lu, Number of WF : %d ..., adcCut : %.0lf, ntracks : %d", 100.0*hipo_nEvent/hipo_nEventMax, hipo_nEvent+1, hipo_nEventMax, nWF, adcCut, trackBank.getRows()).Data());
		return true;
	}
	else {return false;}
}

void Window::drawWaveforms() {
	Grid_waveforms.remove_column(2);
        Grid_waveforms.remove_column(1);
	int num = 0;
	int row = 0;
	int col = 0;
	std::vector<double> bins;
	for (int i = 0; i < NumberOfBins; i++) {
		bins.push_back(1.0*i);
	}
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					double adcMax            = wire->pulse.get_adcMax();
					double adcOffset         = wire->pulse.get_adcOffset();
					double leadingEdgeTime   = wire->pulse.get_leadingEdgeTime();
					double timeOverThreshold = wire->pulse.get_timeOverThreshold();
					double timeMax           = wire->pulse.get_timeMax();
					bool flag = true;
					if (adcMax > 0) {
						if (num >= 10) {break;}
						if (flag_mask_wires) {
							flag = (adcMax >= adcCut) &&
								(adcOffset >= cut_adcOffset_min) &&
								(adcOffset <= cut_adcOffset_max) &&
								(leadingEdgeTime >= cut_leadingEdgeTime_min) &&
								(leadingEdgeTime <= cut_leadingEdgeTime_max) &&
								(timeOverThreshold >= cut_timeOverThreshold_min) &&
								(timeOverThreshold <= cut_timeOverThreshold_max) &&
								(timeMax >= cut_timeMax_min) &&
								(timeMax <= cut_timeMax_max);
							}
						if (!flag) continue;
						std::vector<double> samples = wire->pulse.get_samples();
						// 0 or 1
						col = (num % 2) + 1;
 						row = (num / 2) + 1;
						char buffer[50];
						sprintf(buffer, "L%d W%d", 10*(sl+1) + (l+1), w+1);
						auto area = Gtk::make_managed<Gtk::DrawingArea>();
						area->set_draw_func([this, wire, buffer] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
											cairo_plot_waveform(cr, width, height, wire, buffer);
									      } );
						Grid_waveforms.attach(*area, col, row);
						num++;
					}
				}
			}
		}
	}
	//nWF = num;
}

void Window::drawWaveformsPerLayer() {
	Grid_waveformsPerLayer.remove_column(4);
	Grid_waveformsPerLayer.remove_column(3);
	Grid_waveformsPerLayer.remove_column(2);
	Grid_waveformsPerLayer.remove_column(1);
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				AhdcLayer* layer = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l);
				std::vector<double> sum_samples(NumberOfBins, 0.0); // sum of signals per layer
				double ymin = 0, ymax = 0;
				int nhits = 0;
				for (int w = 0; w < layer->GetNumberOfWires(); w++){
					AhdcWire* wire = layer->GetWire(w);
					double adcMax            = wire->pulse.get_adcMax();
					double adcOffset         = wire->pulse.get_adcOffset();
					double leadingEdgeTime   = wire->pulse.get_leadingEdgeTime();
					double timeOverThreshold = wire->pulse.get_timeOverThreshold();
					double timeMax           = wire->pulse.get_timeMax();
					bool flag = true;
					if (adcMax > 0) {
						if (flag_mask_wires) {
							flag = (adcMax >= adcCut) &&
								(adcOffset >= cut_adcOffset_min) &&
								(adcOffset <= cut_adcOffset_max) &&
								(leadingEdgeTime >= cut_leadingEdgeTime_min) &&
								(leadingEdgeTime <= cut_leadingEdgeTime_max) &&
								(timeOverThreshold >= cut_timeOverThreshold_min) &&
								(timeOverThreshold <= cut_timeOverThreshold_max) &&
								(timeMax >= cut_timeMax_min) &&
								(timeMax <= cut_timeMax_max);
							}
						if (!flag) continue;
						nhits++;
						std::vector<double> samples = wire->pulse.get_samples();
						for (int i = 0; i < NumberOfBins; i++) {
							double value = samples[i];
							sum_samples[i] = sum_samples[i] + value;
							ymin = (ymin < value) ? ymin : value;
							ymax = (ymax > value) ? ymax : value;
						}
					}
				}
				// normalisation
				if (nhits > 0) {
					for (int i = 0; i < NumberOfBins; i++) {
						sum_samples[i] = sum_samples[i]/nhits;
					}
				}
				// Drawings
				auto button = Gtk::make_managed<Gtk::Button>();
				auto area = Gtk::make_managed<Gtk::DrawingArea>();
				auto draw_function = [this, sl, l, layer, ymin, ymax, sum_samples, nhits] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) { // lambda function
					// Define main canvas
					fCanvas canvas(width, height, 0, samplingTime*(NumberOfBins-1), ymin, ymax);
					canvas.define_coord_system(cr);
					canvas.draw_title(cr, "");
					canvas.draw_xtitle(cr, "time (ns)");
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
					// Draw each signals
					int nth = 0;
					for (int w = 0; w < layer->GetNumberOfWires(); w++){
						AhdcWire* wire = layer->GetWire(w);
						double adcMax            = wire->pulse.get_adcMax();
						double adcOffset         = wire->pulse.get_adcOffset();
						double leadingEdgeTime   = wire->pulse.get_leadingEdgeTime();
						double timeOverThreshold = wire->pulse.get_timeOverThreshold();
						double timeMax           = wire->pulse.get_timeMax();
						bool flag = true;
						if (adcMax > 0) {
							if (flag_mask_wires) {
								flag = (adcMax >= adcCut) &&
									(adcOffset >= cut_adcOffset_min) &&
									(adcOffset <= cut_adcOffset_max) &&
									(leadingEdgeTime >= cut_leadingEdgeTime_min) &&
									(leadingEdgeTime <= cut_leadingEdgeTime_max) &&
									(timeOverThreshold >= cut_timeOverThreshold_min) &&
									(timeOverThreshold <= cut_timeOverThreshold_max) &&
									(timeMax >= cut_timeMax_min) &&
									(timeMax <= cut_timeMax_max);
								}
							if (!flag) continue;
							std::vector<double> samples = wire->pulse.get_samples();
							double r1 = 0.0, g1 = 0.0, b1 = 1.0;
							double r2 = 0.0, g2 = 1.0, b2 = 0.0;
							double t = (1.0*nth)/nhits;
							double r = (1-t)*r1 + t*r2;
							double g = (1-t)*g1 + t*g2;
							double b = (1-t)*b1 + t*b2;
							cr->set_source_rgb(r, g, b);
							cr->set_line_width(0.005*seff);
							cr->move_to(x2w(0),y2h(samples[0]));
							for (int i = 1; i < NumberOfBins; i++) {
								// draw a line between points i and i-1
								cr->line_to(x2w(samplingTime*i),y2h(samples[i]));
							}
							cr->stroke();
							nth++;
						}
					}
					// draw the sum of signals
					cr->set_source_rgb(1.0, 0.0, 0.0);
					cr->set_line_width(0.005*seff);
					cr->move_to(x2w(0),y2h(sum_samples[0]));
					for (int i = 1; i < NumberOfBins; i++) {
						// draw a line between points i and i-1
						cr->line_to(x2w(samplingTime*i),y2h(sum_samples[i]));
					}
					cr->stroke();
					// draw frame and axis
					canvas.set_frame_line_width(0.005);
					canvas.draw_frame(cr);
					// add layer name
					char buffer[50];
					sprintf(buffer, "Layer %d", 10*(sl+1) + (l+1));
					cr->set_source_rgb(1.0, 0.0, 0.0);
					cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
					cr->set_font_size(0.6*canvas.get_top_margin());
					Cairo::TextExtents te;
					cr->get_text_extents(buffer, te);
					cr->move_to(0.5*weff - 0.5*te.width, -heff-0.2*canvas.get_top_margin());
					cr->show_text(buffer);
				};
				area->set_draw_func(draw_function);
				button->set_child(*area);
				button->signal_clicked().connect([this, draw_function] () -> void {
							auto window = Gtk::make_managed<Gtk::Window>();
							window->set_title("");
							window->set_default_size(1200,800);
							auto area_bis = Gtk::make_managed<Gtk::DrawingArea>();
							area_bis->set_draw_func(draw_function);
							window->set_child(*area_bis);
							window->show();
						});
				int num_layer = 10*(sl+1) + (l+1);
				int index = -1;
				switch (num_layer) {
					case 11 :
						index = 0;
						break;
					case 21 :
						index = 1;
						break;
					case 22 : 
						index = 2;
						break;
					case 31 : 
						index = 3;
						break;
					case 32 :
						index = 4;
						break;
					case 41 :
						index = 5;
						break;
					case 42 :
						index = 6;
						break;
					case 51 : 
						index = 7;
						break;
					default :
						index = -1;
				}
				int row = (index / 4) + 1 ;
				int col = (index % 4) + 1;
				Grid_waveformsPerLayer.attach(*button, col, row);
			}
		}
	}	
	// end
	Grid_waveformsPerLayer.queue_draw();
}

void Window::drawHistograms() {
	Grid_histograms.remove_column(3);
	Grid_histograms.remove_column(2);
	Grid_histograms.remove_column(1);
	// area 1 : hist1d_adcMax
	auto button1 = Gtk::make_managed<Gtk::Button>();
	auto area1 = Gtk::make_managed<Gtk::DrawingArea>();
	button1->set_child(*area1);
	auto draw_function1 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_adcMax.set_fill_color({0.251, 1, 0.788}); // green
		this->hist1d_adcMax.draw_with_cairo(cr, width, height);
	};
	area1->set_draw_func(draw_function1);
	button1->signal_clicked().connect([this, draw_function1] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area1_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area1_bis->set_draw_func(draw_function1);
				window->set_child(*area1_bis);
				window->show();
			});
	Grid_histograms.attach(*button1,1,1);

	// area 5 : hist1d_adcOffset
	auto button5 = Gtk::make_managed<Gtk::Button>();
	auto area5 = Gtk::make_managed<Gtk::DrawingArea>();
	button5->set_child(*area5);
	auto draw_function5 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_adcOffset.set_fill_color({0.969, 0.78, 0.494}); // orange
		this->hist1d_adcOffset.draw_with_cairo(cr, width, height);
	};
	area5->set_draw_func(draw_function5);
	button5->signal_clicked().connect([this, draw_function5] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area5_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area5_bis->set_draw_func(draw_function5);
				window->set_child(*area5_bis);
				window->show();
			});
	Grid_histograms.attach(*button5,2,1);
	// area 2 : hist1d_leadindEdgeTime
	auto button2 = Gtk::make_managed<Gtk::Button>();
	auto area2 = Gtk::make_managed<Gtk::DrawingArea>();
	button2->set_child(*area2);
	auto draw_function2 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_leadingEdgeTime.set_fill_color({0.961, 0.953, 0.608}); // yellow
		this->hist1d_leadingEdgeTime.draw_with_cairo(cr, width, height);
	};
	area2->set_draw_func(draw_function2);
	button2->signal_clicked().connect([this, draw_function2] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area2_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area2_bis->set_draw_func(draw_function2);
				window->set_child(*area2_bis);
				window->show();
			});
	Grid_histograms.attach(*button2,3,1);
	// area 3 : hist1d_timeOverThreshold
	auto button3 = Gtk::make_managed<Gtk::Button>();
	auto area3 = Gtk::make_managed<Gtk::DrawingArea>();
	button3->set_child(*area3);
	auto draw_function3 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_timeOverThreshold.set_fill_color({0.922, 0.435, 0.647}); // pink (rose)
		this->hist1d_timeOverThreshold.draw_with_cairo(cr, width, height);
	};
	area3->set_draw_func(draw_function3);
	button3->signal_clicked().connect([this, draw_function3] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area3_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area3_bis->set_draw_func(draw_function3);
				window->set_child(*area3_bis);
				window->show();
			});
	Grid_histograms.attach(*button3,1,2);
	// area 4 : hist1d_timeMax
	auto button4 = Gtk::make_managed<Gtk::Button>();
	auto area4 = Gtk::make_managed<Gtk::DrawingArea>();
	button4->set_child(*area4);
	auto draw_function4 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_timeMax.set_fill_color({0.431, 0.765, 0.922}); // blue
		this->hist1d_timeMax.draw_with_cairo(cr, width, height);
	};
	area4->set_draw_func(draw_function4);
	button4->signal_clicked().connect([this, draw_function4] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area4_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area4_bis->set_draw_func(draw_function4);
				window->set_child(*area4_bis);
				window->show();
			});
	Grid_histograms.attach(*button4,2,2);
	// area 6 : hist1d_constantFractionTime
	auto button6 = Gtk::make_managed<Gtk::Button>();
	auto area6 = Gtk::make_managed<Gtk::DrawingArea>();
	button6->set_child(*area6);
	auto draw_function6 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		this->hist1d_constantFractionTime.set_fill_color({0.855, 0.6, 0.969}); // violet
		this->hist1d_constantFractionTime.draw_with_cairo(cr, width, height);
	};
	area6->set_draw_func(draw_function6);
	button6->signal_clicked().connect([this, draw_function6] () -> void {
				auto window = Gtk::make_managed<Gtk::Window>();
				window->set_title("");
				window->set_default_size(1200,800);
				auto area6_bis = Gtk::make_managed<Gtk::DrawingArea>();
				area6_bis->set_draw_func(draw_function6);
				window->set_child(*area6_bis);
				window->show();
			});
	Grid_histograms.attach(*button6,3,2);
	Grid_histograms.queue_draw();
}

void Window::on_mouse_clicked (int n_press, double x, double y) {
	// ************ This portion of code must correspond to the one of on_draw_event !!! ***************
	int width = DrawingArea_event.get_width();
	int height = DrawingArea_event.get_height();
	fCanvas canvas(width, height, -80, 80, -80, 80);
	int window_size = std::min(width,height);
	canvas.set_top_margin(0.05*window_size);
	canvas.set_bottom_margin(0.05*window_size);
	canvas.set_left_margin(0.05*window_size);
	//canvas.set_right_margin(0.05*window_size);
	canvas.set_right_margin(0.15*window_size);
	canvas.set_frame_line_width(0.005);
	canvas.set_stick_width(0.002);
	canvas.set_label_size(0.5);
	canvas.set_title_size(0.6);
	canvas.set_x_start(-80);
	canvas.set_x_end(80);
	canvas.set_y_start(-80);
	canvas.set_y_end(80);
	// **********************             end of the portion of code     *******************************
	auto x2w = [canvas] (double x) {
		return canvas.x2w(x);
	};
	auto y2h = [canvas] (double y) {
		return canvas.y2h(y);
	};
	auto w2x = [canvas] (double w) {
		return canvas.w2x(w);
	};
	auto h2y = [canvas] (double h) {
		return canvas.h2y(h);
	};	
	//printf("n_press : %d, x : %lf, y : %lf\n", n_press, x, y); // x and y coord in widget allocation coordinates
	double w = x - canvas.get_left_margin(); // due to the effect of fCanvas::define_coord_system()
	double h = y - canvas.get_top_margin() - canvas.get_heff(); // due to the effect of fCanvas::define_coord_system()
	double ahdc_x = w2x(w);
	double ahdc_y = h2y(h);	
	//printf("              ahdc_x : %lf, ahdc_y : %lf\n", ahdc_x, ahdc_y);
	// Now, I need to retrieve the closest AhdcWire arond this point (ahdc_x, ahdc_y)
	// then retrieve the layer id
	// then retrieve the waveform
	int layer = -1;
	int component = -1;
	AhdcWire* wire = this->getNearestAhdcWire(ahdc_x, ahdc_y, layer, component);
	if (wire != nullptr) {
		printf("======     AhdcWire  ======\n");
		printf("    x     : %lf\n", wire->x);
		printf("    y     : %lf\n", wire->y);
		printf("    layer : %d\n", layer);
		printf("    comp  : %d\n", component);
		int crate, slot, channel, hv, sub_hv;
		Get_HV_sector(1, layer, component, crate, slot, channel, hv, sub_hv);
		printf("    crate : %d\n", crate);
		printf("    slot  : %d\n", slot);
		printf("    conn  : %d\n", channel/64 + 1);
		printf("    chan  : %d\n", channel);
		printf("    HV_SECTOR : %d-%d\n", hv, sub_hv);
		int tab_number = Book.get_current_page();
		if (tab_number == 0) {
			// popup window
			auto window = Gtk::make_managed<Gtk::Window>();
			window->set_title("AHDC pulse");
			window->set_default_size(1200,800);
			char buffer[50];
			sprintf(buffer, "L%d W%d", layer, component);
			auto area = Gtk::make_managed<Gtk::DrawingArea>();
			area->set_draw_func([this, wire, buffer] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
								cairo_plot_waveform(cr, width, height, wire, buffer);
						      } );
			window->set_child(*area);
			window->show();
		}
		/*char buffer[100];
		sprintf(buffer, "../ressources/get_hv %d %d %d", layer/10, layer % 10, component);
		system(buffer);*/
		TextBuffer_occupancy->set_text("");
		auto iter = TextBuffer_occupancy->end();
		iter = TextBuffer_occupancy->insert(iter, "-----------------------------------\n");
		iter = TextBuffer_occupancy->insert(iter, "|             Output               |\n");
		iter = TextBuffer_occupancy->insert(iter, "-----------------------------------\n");
		iter = TextBuffer_occupancy->insert(iter, "\n");
		{ // x
			char buffer[50];
			sprintf(buffer, "x        :   %lf\n", wire->x);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // y
			char buffer[50];
			sprintf(buffer, "y        :   %lf\n", wire->y);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // layer
			char buffer[50];
			sprintf(buffer, "layer    :   %d\n", layer);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // comp
			char buffer[50];
			sprintf(buffer, "comp     :   %d\n", component);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // crate
			char buffer[50];
			sprintf(buffer, "crate    :   %d\n", crate);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // slot
			char buffer[50];
			sprintf(buffer, "slot     :   %d\n", slot);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // connector
			char buffer[50];
			sprintf(buffer, "connector:   %d\n", channel/64 + 1);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // chan
			char buffer[50];
			sprintf(buffer, "channel   :   %d\n", channel);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		{ // HV sector
			char buffer[50];
			sprintf(buffer, "HV Sector :   %d-%d\n", hv, sub_hv);
			iter = TextBuffer_occupancy->insert(iter, buffer);
		}
		TextView_occupancy.set_buffer(TextBuffer_occupancy);
	}
}

AhdcWire* Window::getNearestAhdcWire(double x, double y, int & _layer, int & _component) {
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					double dx = wire->x - x;
					double dy = wire->y - y;
					double dr = sqrt(dx*dx + dy*dy);
					if (dr < 2) {
						_layer = 10*(sl+1) + (l+1);
						_component = w+1;
					       return wire;	
					} 
				}
			}
		}
	}
	return nullptr;
}

void Window::clearAhdcData() {
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					wire->pulse.reset();
					if (is_reset) wire->occ = 0;
				}
			}
		}
	}
}

void Window::getStats(double & MIN_ADC, double & MAX_ADC, int & MIN_OCC, int & MAX_OCC) {
	AhdcWire* wire0 = ahdc->GetSector(0)->GetSuperLayer(0)->GetLayer(0)->GetWire(0);
	double adc_min = wire0->pulse.get_adcMax();
	double adc_max = adc_min;
	int occ_min = wire0->occ;
	int occ_max = occ_min;
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					double adc = wire->pulse.get_adcMax();
					double occ = wire->occ;
					if (adc < adc_min) adc_min = adc;
					if (adc > adc_max) adc_max = adc;
					if (occ < occ_min) occ_min = occ;
					if (occ > occ_max) occ_max = occ;
				}
			}
		}
	}
	MIN_ADC = adc_min;
	MAX_ADC = adc_max;
	MIN_OCC = occ_min;
	MAX_OCC = occ_max;
}

bool Window::is_oscillating(std::vector<double> samples) {
	int Npts = samples.size();
	if (Npts < 1){ return false; }
	// compute the derived function : df/dx
	std::vector<double> dsamples(Npts, 0.0);
	for (int i = 1; i < Npts; i++) {
		double slope = samples[i] - samples[i-1]; // devided by 1
		dsamples[i-1] = slope;
	}
	dsamples[Npts-1] = 0.0;
	// find the number of zeros
	int nzeros = 0;
	bool flag_skip_next_itr = false;
	for (int i = 2; i < Npts-1; i++) { // bin 0, 1 and Npts-1 are ignored
		if (flag_skip_next_itr) { 
			flag_skip_next_itr = false; 
			continue;
		}
		if (dsamples[i]*dsamples[i-1] < 0) {
			if (dsamples[i+1]*dsamples[i-2] < 0) { // strong condition
				nzeros++;
			}
		}
		else if (dsamples[i]*dsamples[i-1] == 0) {
			if (dsamples[i+1]*dsamples[i-2] < 0) { // strong condition
				nzeros++;
			}
			flag_skip_next_itr = true; // useful to do not count the same zero two times 
		}
		else {
			// do nothing
		}
	}
	//printf("nzeros : %d\n", nzeros);
	// output
	return (nzeros  <= 1);
}
void Window::on_zpos_value_changed() {
	const double val = this->Adjustment_zpos->get_value();
	this->zpos = val;
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
					wire->set_z(this->zpos);
				}
			}
		}
	}
	//printf("AHDC is now view in the plane z = %lf\n", zpos);
	DrawingArea_event.queue_draw();
	drawOccupancy();
}

void Window::drawOccupancy() {
	Grid_occupancy.remove_column(2);
	Grid_occupancy.remove_column(1);
	auto area1 = Gtk::make_managed<Gtk::DrawingArea>();
	area1->set_draw_func(sigc::mem_fun(*this, &Window::on_draw_occupancy) );
	Grid_occupancy.attach(*area1, 1, 1, 1, 2);
	// ...
	auto mouse_click = Gtk::GestureClick::create();
	mouse_click->signal_pressed().connect(sigc::mem_fun(*this, &Window::on_mouse_clicked));
	area1->add_controller(mouse_click);
	// hsitogram 2D
	auto area2 = Gtk::make_managed<Gtk::DrawingArea>();
	auto draw_function2 = [this] (const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
		hist2d_occupancy.set_xtitle("wire");
		hist2d_occupancy.set_ytitle("layer");
		this->hist2d_occupancy.draw_with_cairo(cr, width, height);
	};
	area2->set_draw_func(draw_function2);
	Grid_occupancy.attach(*area2, 2, 1, 1, 1); // left, top, width, height (width and height are used to defined the number during a fusion)
	// Contral
	auto VBox_control = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::VERTICAL, 10);
	Grid_occupancy.attach(*VBox_control, 2,2);
	auto HBox_sectors = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL);
	HBox_sectors->set_margin(10);
	HBox_sectors->set_margin_start(20);
	VBox_control->append(*HBox_sectors);
	// Clear
	auto Button_clear = Gtk::make_managed<Gtk::Button>("Clear");
	HBox_sectors->append(*Button_clear);
	Button_clear->set_hexpand();
	Button_clear->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = -1;
		drawOccupancy();
		printf("Clear HV SECTOR hightlight ...\n");
	});
	// S1
	auto Button_S1 = Gtk::make_managed<Gtk::Button>("S1");
	HBox_sectors->append(*Button_S1);
	Button_S1->set_hexpand();
	Button_S1->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 1;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S2
	auto Button_S2 = Gtk::make_managed<Gtk::Button>("S2");
	HBox_sectors->append(*Button_S2);
	Button_S2->set_hexpand();
	Button_S2->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 2;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S3
	auto Button_S3 = Gtk::make_managed<Gtk::Button>("S3");
	HBox_sectors->append(*Button_S3);
	Button_S3->set_hexpand();
	Button_S3->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 3;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S4
	auto Button_S4 = Gtk::make_managed<Gtk::Button>("S4");
	HBox_sectors->append(*Button_S4);
	Button_S4->set_hexpand();
	Button_S4->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 4;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S5
	auto Button_S5 = Gtk::make_managed<Gtk::Button>("S5");
	HBox_sectors->append(*Button_S5);
	Button_S5->set_hexpand();
	Button_S5->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 5;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S6
	auto Button_S6 = Gtk::make_managed<Gtk::Button>("S6");
	HBox_sectors->append(*Button_S6);
	Button_S6->set_hexpand();
	Button_S6->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 6;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S7
	auto Button_S7 = Gtk::make_managed<Gtk::Button>("S7");
	HBox_sectors->append(*Button_S7);
	Button_S7->set_hexpand();
	Button_S7->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 7;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S8
	auto Button_S8 = Gtk::make_managed<Gtk::Button>("S8");
	HBox_sectors->append(*Button_S8);
	Button_S8->set_hexpand();
	Button_S8->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 8;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// S9
	auto Button_S9 = Gtk::make_managed<Gtk::Button>("S9");
	HBox_sectors->append(*Button_S9);
	Button_S9->set_hexpand();
	Button_S9->signal_clicked().connect([this] () -> void {
		this->HV_SECTOR = 9;
		drawOccupancy();
		printf("Hightlight HV_SECTOR %d ...\n", this->HV_SECTOR);
	});
	// TextView
	VBox_control->append(TextView_occupancy);
	TextView_occupancy.set_margin_start(20);
	TextBuffer_occupancy->set_text("\n");
	auto iter = TextBuffer_occupancy->end();
	iter = TextBuffer_occupancy->insert(iter, "-----------------------------------\n");
	iter = TextBuffer_occupancy->insert(iter, "|             Output               |\n");
	iter = TextBuffer_occupancy->insert(iter, "-----------------------------------\n");
	iter = TextBuffer_occupancy->insert(iter, "\n");
	TextView_occupancy.set_buffer(TextBuffer_occupancy);
}

void Window::on_draw_occupancy(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	fCanvas canvas(width, height, -80, 80, -80, 80);
	int window_size = std::min(width,height);
	canvas.set_top_margin(0.05*window_size);
	canvas.set_bottom_margin(0.05*window_size);
	canvas.set_left_margin(0.05*window_size);
	//canvas.set_right_margin(0.05*window_size);
	canvas.set_right_margin(0.15*window_size);
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
	char face_name[50];
	sprintf(face_name, "(x,y) view in plane z = %.2lf mm", zpos);
	canvas.draw_title(cr, face_name);
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
	// Draw axis orientation
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->set_line_width(0.003*seff);
		// x axis
	cr->move_to(x2w(-70), y2h(-70));
	cr->line_to(x2w(-70+20), y2h(-70));
	cr->stroke();
	cr->move_to(x2w(-70+20), y2h(-70));
	cr->line_to(x2w(-70+20) + 20*cos(M_PI-0.25*M_PI), y2h(-70) + 10*sin(M_PI-0.25*M_PI));
	cr->line_to(x2w(-70+20) + 20*cos(M_PI+0.25*M_PI), y2h(-70) + 10*sin(M_PI+0.25*M_PI));
	cr->close_path();
	cr->fill();
	{ // x axis label
		cr->set_source_rgb(0.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.2*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("x", te);
		cr->move_to(x2w(-70+20)+3, y2h(-70)+3);
		cr->show_text("x");
	}
		// y axis
	cr->move_to(x2w(-70), y2h(-70));
	cr->line_to(x2w(-70), y2h(-70+20));
	cr->stroke();
	cr->move_to(x2w(-70), y2h(-70+20));
	cr->line_to(x2w(-70) + 10*cos(M_PI+0.25*M_PI), y2h(-70+20) - 20*sin(M_PI+0.25*M_PI)); // minus because of the orientation :-(
	cr->line_to(x2w(-70) + 10*cos(-0.25*M_PI), y2h(-70+20) - 20*sin(-0.25*M_PI));
	cr->close_path();
	cr->fill();
	{ // y axis label
		cr->set_source_rgb(0.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.2*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("y", te);
		cr->move_to(x2w(-70)-3-te.width, y2h(-70+20)+3);
		cr->show_text("y");
	}
		// z axis
	cr->set_source_rgb(1.0, 0.0, 0.0);
	cr->move_to(x2w(-70)+10, y2h(-70));
	cr->arc(x2w(-70),y2h(-70), 10, 0, 2*M_PI);
	cr->stroke();
	cr->move_to(x2w(-70)+3, y2h(-70));
	cr->arc(x2w(-70),y2h(-70), 3, 0, 2*M_PI);
	cr->fill();
	{ // z axis label
		cr->set_source_rgb(1.0, 0.0, 0.0);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(1.5*canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents("z", te);
		cr->move_to(x2w(-70) + 13, y2h(-70) + 13 + te.height);
		cr->show_text("z");
	}
	
	// Draw color palette (zbar)
	double adc_min, adc_max;
	int i_occ_min, i_occ_max;
	this->getStats(adc_min, adc_max, i_occ_min, i_occ_max);
	double occ_min = (100.0*i_occ_min)/hipo_nEvent; // expressed in %
	double occ_max = (100.0*i_occ_max)/hipo_nEvent;
	//printf("adc_min : %lf, adc_max : %lf, occ_min : %d, occ_max : %d\n", adc_min, adc_max, occ_min, occ_max);
	fColorPalette Palette(5, 2);
	cr->set_source_rgb(0.0, 0.0, 0.0);
	cr->set_line_width(0.005*seff);
	cr->move_to(x2w(80) + 0.01*window_size, y2h(-80)); // 0.1*window_size is the extra margin of right_margin
	cr->line_to(x2w(80) + 0.1*window_size, y2h(-80));
	cr->line_to(x2w(80) + 0.1*window_size, y2h(80));
	cr->stroke();
	for (int i = 0; i < Palette.get_ncolors(); i++) {
		int ncolors = Palette.get_ncolors();
		fColor rgb = Palette.get_color(i);
		double slope = 160.0/(ncolors); // == 80 - (-80)
		double z = slope*(i-0.0) - 80.0;
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->move_to(x2w(80) + 0.01*window_size, y2h(z));
		cr->line_to(x2w(80) +  0.1*window_size, y2h(z));
		cr->line_to(x2w(80) +  0.1*window_size, y2h(z+1.0*slope));
		cr->line_to(x2w(80) + 0.01*window_size, y2h(z+1.0*slope));
		cr->close_path();
		cr->fill();
	}
	{ // Draw occ_max
		char buffer[50];
		sprintf(buffer, "%.0lf %%", occ_max);
		fColor rgb = Palette.get_color(Palette.get_ncolors()-1);
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents(buffer, te);
		cr->move_to(x2w(80) + 0.05*window_size - 0.5*te.width, y2h(80) - 0.1*canvas.get_top_margin());
		cr->show_text(buffer);
	}
	// Draw occ_min
	{ 
		char buffer[50];
		sprintf(buffer, "%.0lf %%", occ_min);
		fColor rgb = Palette.get_color(0);
		cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
		cr->select_font_face("@cairo:sans-serif",Cairo::ToyFontFace::Slant::NORMAL,Cairo::ToyFontFace::Weight::NORMAL);
		cr->set_font_size(canvas.get_label_size());
		Cairo::TextExtents te;
		cr->get_text_extents(buffer, te);
		cr->move_to(x2w(80) + 0.05*window_size - 0.5*te.width, y2h(-80) + 0.1*canvas.get_top_margin() + te.height);
		cr->show_text(buffer);
	}	
	// Draw AHDC geometry	
	for (int s = 0; s < ahdc->GetNumberOfSectors(); s++) {
		for (int sl = 0; sl < ahdc->GetSector(s)->GetNumberOfSuperLayers(); sl++){
			for (int l = 0; l < ahdc->GetSector(s)->GetSuperLayer(sl)->GetNumberOfLayers(); l++){
				for (int w = 0; w < ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetNumberOfWires(); w++){
					AhdcWire* wire = ahdc->GetSector(s)->GetSuperLayer(sl)->GetLayer(l)->GetWire(w);
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
					//printf("%lf %lf %lf\n", wire->x, wire->y, wire->z);
					cr->move_to(x2w(wire->x) + marker_size, y2h(wire->y)); // correspond to angle == 0 
					cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
					cr->stroke();
					//if ((wire->x == 0) && (wire->y < 0)) {
					if (w == 0) {
						// these wires have a component id == 1 (it is the start of the numerotation)
						cr->set_line_width(0.002*seff);
						cr->set_source_rgba(1.0, 0.0, 0.0, 1.0);
						cr->move_to(x2w(wire->x) + marker_size, y2h(wire->y));
						cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
						cr->stroke();
					}
					// show activated wires
					double occ = wire->occ;
					if (occ > 0) {
						if (occ_min != occ_max) {
							int ncolors = Palette.get_ncolors();
							double occ_percent = (100.0*occ)/hipo_nEvent;
							double slope = (1.0*(ncolors-1))/(occ_max - 0.0);
							int color = floor(slope*(occ_percent-0.0) + 0.0);
							//printf("occ_min : %.0lf, occ_max : %.0lf\n", occ_min, occ_max);
							//printf("color : %d\n", color);
							fColor rgb = Palette.get_color(color);
							cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
						} else {
							fColor rgb = Palette.get_color(0);
							cr->set_source_rgb(rgb.r, rgb.g, rgb.b);
						}
						double marker_size = std::min(2.0*weff/(x_end-x_start), 2.0*heff/(y_end-y_start));
						cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
						cr->fill();
					}
					// highlight HV_SECTOR
					int layer = 10*(sl+1) + (l+1);
					int component = w+1;
					int crate, slot, channel, hv, sub_hv;
					Get_HV_sector(1, layer, component, crate, slot, channel, hv, sub_hv);
					if (HV_SECTOR == hv) {
						if 	(sub_hv == 1) {
							cr->set_source_rgba(127/255.0, 0/255.0, 255/255.0, 1.0);
						}
						else if (sub_hv == 2) {
							cr->set_source_rgba(1.0, 0.0, 0.945, 1.0);
						}
						else if (sub_hv == 3) {
							cr->set_source_rgba( 1.0, 0.0, 0.0, 1.0);
						}
						cr->arc(x2w(wire->x), y2h(wire->y) , marker_size, 0, 2*M_PI);
						cr->fill();
					}
				
				}
			}
		}
	}
	// TEST
	/*TCanvas* canvas1 = new TCanvas("all","all hist 2d",1800, 990);
	hist2d_occupancy.get_root_histo()->Draw("COLZ");
	canvas1->Print("../output/test_fh2d.pdf");
	delete canvas1;	*/
}

void Window::on_draw_3Dview(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height) {
	//printf("w : %d, h : %d\n", width, height);
	//cr->translate(0, height);
	//cr->scale(1.0, -1.0);
	//cr->move_to(50,50);
	//cr->set_source_rgb(0.0, 0.0, 0.0);
	//cr->set_line_width(10);
	//cr->rectangle(80,20, 100, 100);
	//cr->stroke();
	Frame::test1(angle_alpha, angle_beta, angle_gamma, cr, width, height);	
}

int Window::layer2number(int digit) {
	if      (digit == 11) {
		return 1;
	} 
	else if (digit == 21) {
		return 2;
	} 
	else if (digit == 22) {
		return 3;
	} 
	else if (digit == 31) {
		return 4;
	} 
	else if (digit == 32) {
		return 5;
	} 
	else if (digit == 41) {
		return 6;
	} 
	else if (digit == 42) {
		return 7;
	} 
	else if (digit == 51) {
		return 8;
	} else {
		return -1; // not a layer
	}
}

void Window::Get_HV_sector(int sector, int layer, int component, int & crate, int & slot, int & channel, int & hv, int & sub_hv) {
	AhdcMapping::GetDreamChannel(1, layer, component, crate, slot, channel);
	hv = -1;
	sub_hv = -11;
	// read hv_notes.txt
	if (slot == 1) { 
		hv = (channel-1)/64 + 1;
		if (hv % 2 == 0) hv--;
		else hv++;
		int num = (channel-1)%64 + 1;
		if ((num >= 1) && (num <= 23)) {
			sub_hv = 3;
		}
		else if ((num >= 24) && (num <= 45)) {
			sub_hv = 2;
		}
		else if ((num >= 46) && (num <= 64)) {
			sub_hv = 1;
		}
	}
	else if (slot == 4) { // hv only 9
		// read hv_notes.txt
		hv = (channel-1)/64 + 1;
		if (hv == 2) { hv = 9;}
		if (hv != 9) {
			if (hv % 2 == 0) hv--;
			else hv++;
		}
		int num = (channel-1)%64 + 1;
		if ((num >= 1) && (num <= 23)) {
			sub_hv = 3;
		}
		else if ((num >= 24) && (num <= 45)) {
			sub_hv = 2;
		}
		else if ((num >= 46) && (num <= 64)) {
			sub_hv = 1;
		}
	}
}

/** Main function */
int main (int argc, char * argv[]) {
	std::cout << "Start GUi..." << std::endl;

	auto app = Gtk::Application::create("org.gtkmm.example");

	return app->make_window_and_run<Window>(argc, argv);
}



