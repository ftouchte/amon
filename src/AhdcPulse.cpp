/****************************************************
 * AhdcPulse.cpp
 *
 * @author Felix Touchte Codjo
 * @date March 25, 2025
 * *************************************************/

#include "AhdcPulse.h"

AhdcPulse::AhdcPulse() {
	integral = 0;
	adcMax = -1;
	adcOffset = 0;
	timeMax = 0;
	leadingEdgeTime = 0;
	timeOverThreshold = 0;
	constantFractionTime = 0;
}

void AhdcPulse::reset() {
	samples.clear();	
	integral = 0;
	adcMax = 0;
	adcOffset = 0;
	timeMax = 0;
	leadingEdgeTime = 0;
	timeOverThreshold = 0;
	constantFractionTime = 0;
}
void AhdcPulse::set_integral(double _integral) { integral = _integral;}
void AhdcPulse::set_adcMax(double _adcMax) { adcMax = _adcMax;}
void AhdcPulse::set_adcOffset(double _adcOffset) { adcOffset = _adcOffset;}
void AhdcPulse::set_timeMax(double _timeMax) { timeMax = _timeMax;}
void AhdcPulse::set_leadingEdgeTime(double _leadingEdgeTime) { leadingEdgeTime = _leadingEdgeTime;}
void AhdcPulse::set_timeOverThreshold(double _timeOverThreshold) { timeOverThreshold = _timeOverThreshold;}
void AhdcPulse::set_constantFractionTime(double _constantFractionTime) { constantFractionTime = _constantFractionTime;}
void AhdcPulse::set_samples(std::vector<double> _samples) { samples = _samples;}

double AhdcPulse::get_integral() {return integral;}
double AhdcPulse::get_adcMax() {return adcMax;}
double AhdcPulse::get_adcOffset() {return adcOffset;}
double AhdcPulse::get_timeMax() {return timeMax;}
double AhdcPulse::get_leadingEdgeTime() {return leadingEdgeTime;}
double AhdcPulse::get_timeOverThreshold() {return timeOverThreshold;}
double AhdcPulse::get_constantFractionTime() {return constantFractionTime;}
std::vector<double> AhdcPulse::get_samples() {return samples;}

