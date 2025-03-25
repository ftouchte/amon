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
	double integral;
	double adcMax;
	double adcOffset;
	double timeMax;
	double leadingEdgeTime;
	double timeOverThreshold;
	double constantFractionTime;
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
	void set_samples(std::vector<double> _samples);

	double get_integral();
	double get_adcMax();
	double get_adcOffset();
	double get_timeMax();
	double get_leadingEdgeTime();
	double get_timeOverThreshold();
	double get_constantFractionTime();
	std::vector<double> get_samples();
};

#endif
