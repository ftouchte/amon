/*************************************************
 * Class for 3D Point
 *
 * @author Felix Touchte Codjo
 * @date November 18,2024
 * **********************************************/

#ifndef FELIX_POINT_3D_H
#define FELIX_POINT_3D_H

class Point3D {
public :
	double x; ///< x coordinates
	double y; ///< y coordinates
	double z; ///< z coordinates

	Point3D(double _x = 0, double _y = 0, double _z = 0);
	Point3D(const Point3D & pt);
	Point3D & operator=(const Point3D & pt);
	void Print();

	void Translate(double a, double b, double c);
};

#endif
