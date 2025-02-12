/***********************************************
 * Class for 1D histogram
 *
 * designed to be used in gtkmm
 * drawing area
 *
 * @author Felix Touchte Codjo
 * ********************************************/
#include <cstdio>
#include <cmath>

fH1D::fH1D(std::string _title = "", int _nbins, double _xmin, double _xmax) : title(_title), nbins(_nbins), xmin(_xmin), xmax(_xmax) {
	if (xmax < xmin) {
		printf("Histogram parameters are incorrects : xmax < xmin");		
		return ;
	}
	binw = (xmax - xmin)/nbins;
	for (int i = 0; i < nbins; i++) {
		binArray.push_back(xmin + i*binw + 0.5*binw);
		binBuffer.push_back(0.0);
	}
	max_buffer = 0;
	min_buffer = 0;
	underflow = 0;
	overflow = 0;
	nEntries = 0;
	sum = 0;
	sum2 = 0;
}

void fH1D::fill(double x) {
	int bin = getBinNumber(x);
	if (bin == -1) { underflow++;}
	if (bin == -11) { overflow++;}
	nEntries++;
	binBuffer[i] = binBuffer[i] + 1.0;
	// stats
	sumw += 1;
	sum += x;
	sum2 += x*x;
}

/**
 * fill x with weight of w
 * @param x value to fill
 * @param w weight
 *
 * @note return -1 if undeflow and -11 if overflow
 */
void fH1D::fill(double x, double w) {
	int bin = getBinNumber(x);
	if (bin == -1) { underflow++;}
	if (bin == -11) { overflow++;}
	nEntries++;
	binBuffer[i] = binBuffer[i] + w;
	//stats
	sumw += w;
	sum += w*x;
	sum2 += w*x*x;
}
/**
 *
 * @note numerotation starts at 0
 */
int fH1D::getBinNumber(double x) const {
	if (x < xmin) {return -1;}
	if (x > xmax) {return -11;}
	for (int i = 0; i < nbins; i++) {
		xinf = xmin + i*binw;
		xsup = xinf + binw;
		if ((x >= xinf) && (x < xsup)) {
			return i;
		}
	}
}

double fH1D::getBinBuffer(int bin) const {
	if ((bin < 0) || (bin >= nbins)) {
		return 0;
	}
	return binBuffer[bin];
}

double fH1D::getBinValue(int bin) const {
	/*if ((bin < 0) || (bin >= nbins)) {
		return 0;
	}*/
	return binArray[bin];
}

int fH1D::getEntries() const { return nEntries;}
double fH1D::getMean() const { return sum/sumw;}
double fh1D::getStDev() const { return sqrt(sum2/sumw - getMean()*getMean());}
double fh1D::getBinWidth() const {return binw;}
int fh1D::getNumberOfBins() const {return nbins;}
std::vector<double> fh1D::getBinArray() const {return binArray;}
std::vector<double> fh1D::getBinBuffer() const {return binBuffer;}
void fh1D::set_xtitle(std::string name) {xtitle = name;}
void fh1D::set_ytitle(std::string name) {ytitle = name;}

double fH1D::getMax() const {
	int vmax = 0;
	for (int i = 0; i < nbins; i++) {
		vmax = (vmax < binBuffer[i]) ? binBuffer[i] : vmax;
	}
	return vmax;
}

//inclure tout gtkmm
//ou télécharger cairo
//draw(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height){
void fH1D::draw() {
	// prévoire les offset
	ax = fAxis(xmin, xmax, 10, 5, 0);
	ay = fAxis(0, getMax(), 10, 5, 0);	
}
