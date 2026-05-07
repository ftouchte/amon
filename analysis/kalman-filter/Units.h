/**********************************
 * Standard units
 * 
 * @author Felix Touchte Codjo
 * @date April 02, 2026
 **********************************/

#ifndef F_UNITS_H
#define F_UNITS_H

#include <cmath>

/**
 * @brief Automatic conversion to default units. 
 * 
 * Energy --> GeV
 * 
 * Distance --> cm
 * 
 * Angles --> rad
 * 
 * Time --> ns
 * 
 */
struct Units {

    // Energy (GeV)
    static constexpr double GeV = 1.0;
    static constexpr double MeV = 1e-3 * GeV;
    static constexpr double keV = 1e-6 * GeV;
    static constexpr double  eV = 1e-9 * GeV;

    // Distance (cm)
    static constexpr double cm = 1.0;
    static constexpr double mm = 1e-1 * cm;
    static constexpr double  m = 1e2 * cm;

    // Angles (rad)
    static constexpr double rad = 1.0;
    static constexpr double deg = M_PI * rad / 180;

    // Times
    static constexpr double ns = 1.0;
    static constexpr double ps = 1e-3 * ns;
    static constexpr double µs = 1e3  * ns;
    static constexpr double ms = 1e6  * ns;
    static constexpr double sec = 1e9 * ns;

};



#endif
