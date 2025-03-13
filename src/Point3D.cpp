/*************************************************
 * Class for 3D Point
 *
 * @author Felix Touchte Codjo
 * @date November 18,2024
 * **********************************************/

#include "Point3D.h"
#include <stdio.h>

/** Constructor */
Point3D::Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

/** Copy constructor */
Point3D::Point3D(const Point3D & pt) {
	x = pt.x;
	y = pt.y;
	z = pt.z;
}

/** Printing */
void Point3D::Print(){
	printf("(%.1lf,%.1lf,%.1lf)\n",x,y,z);
}

/** Apply a linear translation */
void Point3D::Translate(double a, double  b, double c){
	x += a;
	y += b;
	z += c;
}

/** Operator= */
Point3D & Point3D::operator=(const Point3D & pt){
	x = pt.x;
	y = pt.y;
	z = pt.z;
	return *this;
}
