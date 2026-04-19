package io.github.ftouchte.alignment;

import java.util.ArrayList;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;

/**
 * Contains the list of all histos used during the analysis.
 * 
 * @note add new histos (H1F or H2F) and complete the constructor and the merge methods (for concurrent programming)
 */
public class Histos {

    /** 1D residuals LR per layers */
    ArrayList<H1F> h1_residual_LR_per_layers = new ArrayList<>();

    /** 2D residuals LR per layers versus vz*/
    ArrayList<H2F> h2_corr_vz_residual_LR_per_layers = new ArrayList<>();

    /** 1D residuals per layers */
    ArrayList<H1F> h1_residual_per_layers = new ArrayList<>();
    
    /** 1D residuals LR per wire */
    ArrayList<H1F> h1_residual_LR_per_wires = new ArrayList<>();

    /** 2D residuals LR per wires versus vz*/
    ArrayList<H2F> h2_corr_vz_residual_LR_per_wires = new ArrayList<>();

    /** track theta */
    H1F h1_track_theta;

    /** track delta phi with electron*/
    H1F h1_track_delta_phi;

    /** Time2Distance per layer */
    ArrayList<H2F> h2_time2distance_per_layers = new ArrayList<>();

    /** Time2Distance per wires */
    ArrayList<H2F> h2_time2distance_per_wires = new ArrayList<>();


    /**
     * Create an instance of {@link Histos}
     * @param niter Iteration number
     */
    public Histos(int niter) {

        /// h1 residual LR per layers
        for (int i = 0; i < 9; i++) {
            H1F h = new H1F("residual-LR-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, "residual LR (layer " + AhdcWireId.number2layer(i) + ")", 100, -1.5, 1.5);
            h.setTitleX("layer " + AhdcWireId.number2layer(i) + ", residual LR (mm)");
            if (i == 0) h.setTitleX("all layers, residual LR (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_LR_per_layers.add(h);
        }

        /// h1 residual per layers
        for (int i = 0; i < 9; i++) {
            H1F h = new H1F("residual-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, "residual (layer " + AhdcWireId.number2layer(i) + ")", 100, -1.5, 1.5);
            h.setTitleX("residual (mm)");
            if (i == 0) h.setTitleX("all layers, residual (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_per_layers.add(h);
        }

        /// h1 residual LR per wires
        for (int i = 0; i < 576; i++) {
            H1F h = new H1F("residual-LR-wire-" + i + "itr-" + niter, "residual-LR-wire-" + i + "itr-" + niter, 100, -3, 3);
            AhdcWireId id = new AhdcWireId(i);
            h.setTitleX("wire " + i + ", L" + id.layer + "C" + id.component + ", residual LR (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_LR_per_wires.add(h);
        }

        /// h2 corr vz residual LR per layers
        for (int i = 0; i < 9; i++) {
            H2F h = new H2F("corr-vz-residual-LR-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, 60, -215, 160, 120, -1.5, 1.5);
            h.setTitleX("layer " + AhdcWireId.number2layer(i) + ", vz (mm)");
            if (i == 0) h.setTitleX("all layers, vz(mm)");
            if (i % 3 == 0) h.setTitleY("residual LR (mm)");
            h2_corr_vz_residual_LR_per_layers.add(h);
        }

        /// h2 corr vz residual LR per wires
        for (int i = 0; i < 576; i++) {
            H2F h = new H2F("corr-vz-residual-LR-wire-" + i + "itr-" + niter, 60, -215, 160, 120, -1.5, 1.5);
            AhdcWireId identifier = new AhdcWireId(i);
            h.setTitleX("L" + identifier.layer + "W" + identifier.component +  ", vz (mm)");
            if (i % 3 == 0) h.setTitleY("residual LR (mm)");
            h2_corr_vz_residual_LR_per_wires.add(h);
        }

        /// h2 time2distance per layers
        for (int i = 0; i < 9; i++) {
            H2F h = new H2F("time2distance-" + AhdcWireId.number2layer(i) + "itr-" + niter, 100, 0, 250, 100, 0, 3);
            h.setTitleX("layer " + AhdcWireId.number2layer(i) + ", time (ns)");
            if (i == 0) h.setTitleX("all layers, time (ns)");
            if (i % 3 == 0) h.setTitleY("distance (mm)");
            h2_time2distance_per_layers.add(h);
        }

        /// h2 time2distance per wires
        for (int i = 0; i < 576; i++) {
            H2F h = new H2F("time2distance-" + i + "itr-" + niter, 100, 0, 250, 100, 0, 3);
            AhdcWireId identifier = new AhdcWireId(i);
            h.setTitleX("L" + identifier.layer + "W" + identifier.component +  ", time (ns)");
            if (i % 3 == 0) h.setTitleY("distance (mm)");
            h2_time2distance_per_wires.add(h);
        }

        /// h1 track theta
        h1_track_theta = new H1F("track-theta-itr-" + niter, "track theta", 100, 0, 180);
        h1_track_theta.setTitleX("theta (deg)");
        h1_track_theta.setTitleY("count");

        /// h1 track delta phi
        h1_track_delta_phi = new H1F("track-delat-phi-itr-" + niter, "track theta", 100, -10, 10);
        h1_track_delta_phi.setTitleX("delta phi (deg)");
        h1_track_delta_phi.setTitleY("count");
    }

    /** 
     * Add the content of the incoming Histos object int the current one.
     * Very useful for the merge operation after the parallelisation.
     * @param histos
     */
    public void merge(Histos histos) {
        // Per layers
        for (int i = 0; i < 9; i++) {
            // h1 residual LR
            this.h1_residual_LR_per_layers.get(i).add(histos.h1_residual_LR_per_layers.get(i));
            // h1 residual
            this.h1_residual_per_layers.get(i).add(histos.h1_residual_per_layers.get(i));
            // h2 corr vz with residual LR
            this.h2_corr_vz_residual_LR_per_layers.get(i).add(histos.h2_corr_vz_residual_LR_per_layers.get(i));
            // h2 time2distance
            this.h2_time2distance_per_layers.get(i).add(histos.h2_time2distance_per_layers.get(i));
        }
        // Per wires
        for (int i = 0; i < 576; i++) {
            this.h1_residual_LR_per_wires.get(i).add(histos.h1_residual_LR_per_wires.get(i));
            this.h2_corr_vz_residual_LR_per_wires.get(i).add(histos.h2_corr_vz_residual_LR_per_wires.get(i));
            // h2 time2distance
            this.h2_time2distance_per_wires.get(i).add(histos.h2_time2distance_per_wires.get(i));
        }

        // Integrated
        this.h1_track_theta.add(histos.h1_track_theta);
        this.h1_track_delta_phi.add(histos.h1_track_delta_phi);

    }
}
