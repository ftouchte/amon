package io.github.ftouchte.utils;

/**
 * The default unit for the energy is GeV
 * The default unit for the distance is cm
 * The default unit for the angle is rad
 */
public class Units {
    // Energy
    public static final double GeV = 1.0;
    public static final double MeV = 1e-3 * GeV;
    public static final double keV = 1e-6 * GeV;
    public static final double  eV = 1e-9 * GeV;

    // Distance
    public static final double cm = 1.0;
    public static final double mm = 1e-1 * cm;
    public static final double  m = 1e2 * cm;

    // Angles
    public static final double rad = 1.0;
    public static final double deg = Math.PI * rad / 180;

}
