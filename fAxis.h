/************************************************* 
 * Class for axis (labels and sticks)
 * 
 * @author Felix Touchte Codjo 
 * @date January 26,2025 
 * **********************************************/ 
 
#ifndef FELIX_AXIS_H
#define FELIX_AXIS_H

#include <vector>
#include <string>

class fAxis {
private :
	double start; ///< start point
	double end; ///< end point
	double width; ///< length of the interval
	double eps1; ///< epsilon, precision of the first division
	double eps2; ///< precision of the second division
	double eps3; ///< precision of the third division
	int n1; ///< number of first divisions
	int n2; ///< number of second divisions
	int n3; ///< number of third divisions
	std::vector<double> Div1; ///< first divisions
	std::vector<double> Div2; ///< second divisions
	std::vector<double> Div3; ///< third divisions
	std::vector<std::string> Labels1; ///< Labels of the first division (formatted version of Div1)
	std::vector<std::string> Labels2; ///< Labels of the second division (formatted version of Div2)
	std::vector<std::string> Labels3; ///< Labels of the third division (formatted version of Div3)

public :
	fAxis(double _start, double _end, int _n1 = 10, int _n2 = 5, int _n3 = 0);
	~fAxis();
	int get_n1();
	int get_n2();
	int get_n3();
	int get_eps1();
	int get_eps2();
	int get_eps3();
	void print();
	std::vector<std::string> get_labels1();
	std::vector<std::string> get_labels2();
	std::vector<std::string> get_labels3();
};

#endif
