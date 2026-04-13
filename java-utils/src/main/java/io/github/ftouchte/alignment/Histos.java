package io.github.ftouchte.alignment;

import java.util.ArrayList;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;

/**
 * Contains the list of all histos used during the analysis.
 * 
 * @note add new histos (H1F or H2F) and complete the constructor and the merge methods.
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

    /**
     * Create an instance of {@link Histos}
     * @param niter Iteration number
     */
    public Histos(int niter) {
        // h1 residual LR per layers
        for (int i = 0; i < 9; i++) {
            H1F h = new H1F("residual-LR-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, "residual LR (layer " + AhdcWireId.number2layer(i) + ")", 100, -1.5, 1.5);
            h.setTitleX("layer " + AhdcWireId.number2layer(i) + ", residual LR (mm)");
            if (i == 0) h.setTitleX("all layers, residual LR (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_LR_per_layers.add(h);
        }
        // h1 residual
        for (int i = 0; i < 9; i++) {
            H1F h = new H1F("residual-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, "residual (layer " + AhdcWireId.number2layer(i) + ")", 100, -3, 3);
            h.setTitleX("residual (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_per_layers.add(h);
        }
        // h1 residual LR per wires
        for (int i = 0; i < 576; i++) {
            H1F h = new H1F("residual-LR-wire-" + i + "itr-" + niter, "residual-LR-wire-" + i + "itr-" + niter, 100, -3, 3);
            AhdcWireId id = new AhdcWireId(i);
            h.setTitleX("wire " + i + ", L" + id.layer + "C" + id.component + ", residual LR (mm)");
            h.setTitleY("count");
            //h.setOptStat(1111);
            h1_residual_LR_per_wires.add(h);
        }

        /// --- histograms 2D
        for (int i = 0; i < 9; i++) {
            H2F h = new H2F("corr-vz-residual-LR-layer-" + AhdcWireId.number2layer(i) + "itr-" + niter, 60, -16, 16, 120, -1.5, 1.5);
            h.setTitleX("layer " + AhdcWireId.number2layer(i) + ", vz (cm)");
            if (i == 0) h.setTitleX("all layers, vz(cm)");
            if (i % 3 == 0) h.setTitleY("residual LR (mm)");
            h2_corr_vz_residual_LR_per_layers.add(h);
        }
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
        }
        // Per wires
        for (int i = 0; i < 576; i++) {
            this.h1_residual_LR_per_wires.get(i).add(histos.h1_residual_LR_per_wires.get(i));
        }

    }
}
