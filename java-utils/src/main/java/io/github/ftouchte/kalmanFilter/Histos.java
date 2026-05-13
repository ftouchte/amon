package io.github.ftouchte.kalmanFilter;

import org.jlab.groot.data.H1F;

public class Histos {
    /** Reconstructed kinematics */
    H1F h1_p, h1_theta, h1_phi, h1_vz;

    /** Expected kinematics*/
    H1F h1_p0, h1_theta0, h1_phi0, h1_vz0;

    /** Differences between reconstructed and expected distributions */
    H1F h1_delta_p, h1_delta_theta, h1_delta_phi, h1_delta_vz;

    /**
     * Constructor. Place to initialise all histograms
     */
    public Histos() {

        h1_p = new H1F("reconstructed_p", "reconstructed p", 100, 0, 1000);
        h1_theta = new H1F("reconstructed_theta", "reconstructed theta", 100, 0, 180);
        h1_phi = new H1F("reconstructed_phi", "reconstructed phi", 100, 0, 360);
        h1_vz = new H1F("reconstructed_vz", "reconstructed vz", 100, -24, 20);

        h1_p0 = new H1F("expected_p", "expected p", 100, 0, 1000);
        h1_theta0 = new H1F("expected_theta", "expected theta", 100, 0, 180);
        h1_phi0 = new H1F("expected_phi", "expected phi", 100, 0, 360);
        h1_vz0 = new H1F("expected_vz", "expected vz", 100, -24, 20);

        h1_delta_p = new H1F("delta_p", "delta p", 100, 0, 1000);
        h1_delta_theta = new H1F("delta_theta", "delta theta", 100, 0, 180);
        h1_delta_phi = new H1F("delta_phi", "delta phi", 100, 0, 360);
        h1_delta_vz = new H1F("delta_vz", "delta vz", 100, -24, 20);

    }

    /** 
     * Add the content of the incoming Histos object int the current one.
     * Very useful for the merge operation after the parallelisation.
     * @param histos
     */
    public void merge(Histos histos) {

        this.h1_p.add(histos.h1_p);
        this.h1_theta.add(histos.h1_theta);
        this.h1_phi.add(histos.h1_phi);
        this.h1_vz.add(histos.h1_vz);

        this.h1_p0.add(histos.h1_p0);
        this.h1_theta0.add(histos.h1_theta0);
        this.h1_phi0.add(histos.h1_phi0);
        this.h1_vz0.add(histos.h1_vz0);

        this.h1_delta_p.add(histos.h1_delta_p);
        this.h1_delta_theta.add(histos.h1_delta_theta);
        this.h1_delta_phi.add(histos.h1_delta_phi);
        this.h1_delta_vz.add(histos.h1_delta_vz);
        
    }
}
