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
    wfType = 100; // non defined
    mask = true;
    nhits = 0;
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
    wfType = -1; // non defined
    mask = true;
    nhits = 0;
}

void AhdcPulse::set_integral(double _integral) { integral = _integral;}
void AhdcPulse::set_adcMax(double _adcMax) { adcMax = _adcMax;}
void AhdcPulse::set_adcOffset(double _adcOffset) { adcOffset = _adcOffset;}
void AhdcPulse::set_timeMax(double _timeMax) { timeMax = _timeMax;}
void AhdcPulse::set_leadingEdgeTime(double _leadingEdgeTime) { leadingEdgeTime = _leadingEdgeTime;}
void AhdcPulse::set_timeOverThreshold(double _timeOverThreshold) { timeOverThreshold = _timeOverThreshold;}
void AhdcPulse::set_constantFractionTime(double _constantFractionTime) { constantFractionTime = _constantFractionTime;}
void AhdcPulse::set_binOffset(int _binOffset) { binOffset = _binOffset;}
void AhdcPulse::set_samples(std::vector<double> _samples) { samples = _samples;}
void AhdcPulse::set_wfType(int _type) { wfType = _type;}
void AhdcPulse::set_mask(bool _mask) { mask = _mask;}
void AhdcPulse::triggered() { nhits++;}

double AhdcPulse::get_integral() const {return integral;}
double AhdcPulse::get_adcMax() const {return adcMax;}
double AhdcPulse::get_adcOffset() const {return adcOffset;}
double AhdcPulse::get_timeMax() const {return timeMax;}
double AhdcPulse::get_leadingEdgeTime() const {return leadingEdgeTime;}
double AhdcPulse::get_timeOverThreshold() const {return timeOverThreshold;}
double AhdcPulse::get_constantFractionTime() const {return constantFractionTime;}
int    AhdcPulse::get_binOffset() const {return binOffset;}
std::vector<double> AhdcPulse::get_samples() const {return samples;}
int AhdcPulse::get_wfType() const { return wfType;}
bool AhdcPulse::get_mask() const { return mask;}
bool AhdcPulse::is_masked() const { return !mask;}
int AhdcPulse::get_nhits() const { return nhits;}
