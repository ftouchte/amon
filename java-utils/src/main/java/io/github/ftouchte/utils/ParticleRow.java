package io.github.ftouchte.utils;

/**
 * 
 * 
 */
public class ParticleRow {
    private int row = -1; // row in the HIPO bank
    private double px, py, pz; // GeV
    private double p, theta, phi; // GeV, radian, radian
    private double pT;

    private double vx, vy, vz; // cm
    
    /**
     * If no {@link Units#Units} is specified, the unit is assumed to be GeV (energies) and cm (distances)
     * @param _px momentum px component
     * @param _py momentum py component
     * @param _pz momentum pz component
     */
    public ParticleRow(double _px, double _py, double _pz) {
        px = _px;
        py = _py;
        pz = _pz;
        p = Math.sqrt(px*px + py*py + pz*pz);
        pT = Math.sqrt(px*px + py*py);
        theta = Math.acos(pz/p);
        phi = Math.atan2(py, px);
    }

    /** Return the amplitude of the momentum vector */
    public double p() {
        return p;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the momentum accordingly to the given unit
     */
    public double p(double unit) {
        return p/unit;
    }

    /** Return the amplitude of the transverse momentum vector */
    public double pT() {
        return pT;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the of the transverse momentum accordingly to the given unit
     */
    public double pT(double unit) {
        return pT/unit;
    }

    /** Return the theta angle (spherical coordinate) of the momentum vector */
    public double theta() {
        return theta;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the theta angle (spherical coordinate) of the momentum vector component accordingly to the given unit
     */
    public double theta(double unit) {
        return theta/unit;
    }

    /** Return the phi angle (spherical coordinate) of the momentum vector */
    public double phi() {
        return phi;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the theta angle (spherical coordinate) of the momentum vector component accordingly to the given unit
     */
    public double phi(double unit) {
        return phi/unit;
    }
    
    /** Return the px component of the momentum in GeV*/
    public double px() {
        return px;
    }

    /** Return the py component of the momentum in GeV*/
    public double py() {
        return py;
    }

    /** Return the pz component of the momentum in GeV*/
    public double pz() {
        return pz;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the momentum px component accordingly to the given unit
     */
    public double px(double unit) {
        return px/unit;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the momentum py component accordingly to the given unit
     */
    public double py(double unit) {
        return py/unit;
    }

    /**
     * 
     * @param unit should be a {@link Units#Units}
     * @return return the momentum pz component accordingly to the given unit
     */
    public double pz(double unit) {
        return pz/unit;
    }

    /**
     * Set vertex
     * @param _vx
     * @param _vy
     * @param _vz
     */
    public void SetVertex(double _vx, double _vy, double _vz) {
        vx = _vx;
        vy = _vy;
        vz = _vz;
    }

    /** Return vx vertex */
    public double vx() {
        return vx;
    }

    /** Return vy vertex */
    public double vy() {
        return vy;
    }

    /** Return vz vertex */
    public double vz() {
        return vz;
    }

    /**
     * Set HIPO bank row
     * @param _row
     */
    public void SetBankRow(int _row) {
        row = _row;
    }

    /** Return bank row */
    public int GetBankRow() {
        return row;
    }

}
