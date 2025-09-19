/****************************************************
 * AhdcPulse.h
 *
 * @author Felix Touchte Codjo
 * @date March 25, 2025
 * *************************************************/

#ifndef AHDC_PULSE
#define AHDC_PULSE

#include <vector>

class AhdcPulse {
	std::vector<double> samples;
	int wfType;
	double integral;
	double adcMax;
	double adcOffset;
	double timeMax;
	double leadingEdgeTime;
	double timeOverThreshold;
	double constantFractionTime;
	int binOffset;
    bool mask;
    int nhits; // it may happen that we have several hits on this wires, in that case we should have vectors of waveforms, wfType, integral, ... TO BE DONE later
public :
	AhdcPulse();
	void reset();
	void set_integral(double _integral);
	void set_adcMax(double _adcMax);
	void set_adcOffset(double _adcOffset);
	void set_timeMax(double _timeMax);
	void set_leadingEdgeTime(double _leadingEdgeTime);
	void set_timeOverThreshold(double _timeOverThreshold);
	void set_constantFractionTime(double _constantFractionTime);
	void set_binOffset(int _binOffset);
	void set_samples(std::vector<double> _samples);
	void set_wfType(int _type);
    void set_mask(bool _mask);
    void triggered();

	double get_integral() const;
	double get_adcMax() const;
	double get_adcOffset() const;
	double get_timeMax() const;
	double get_leadingEdgeTime() const;
	double get_timeOverThreshold() const;
	double get_constantFractionTime() const;
	int get_binOffset() const;
	std::vector<double> get_samples() const;
    int get_wfType() const;
    bool get_mask() const;
    bool is_masked() const;
    int get_nhits() const;
};

#endif
