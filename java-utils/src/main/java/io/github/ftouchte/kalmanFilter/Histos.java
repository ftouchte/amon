package io.github.ftouchte.kalmanFilter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.jfree.chart.JFreeChart;
import org.jlab.groot.data.H1F;

import com.itextpdf.text.DocumentException;

import io.github.ftouchte.kalmanFilter.Renderer.RendererOutputType;

public class Histos {
    /** Reconstructed kinematics */
    H1F h1_p, h1_theta, h1_phi, h1_vz;

    /** Expected kinematics*/
    H1F h1_p0, h1_theta0, h1_phi0, h1_vz0;

    /** Differences between reconstructed and expected distributions */
    H1F h1_delta_p, h1_delta_theta, h1_delta_phi, h1_delta_vz;

    /** ALERTEngine computing time per event*/
    H1F h1_computing_time;

    /** Residuals */
    H1F h1_residual;

    Map<String, H1F> map_histo1d = new HashMap<>();

    /** Tag */
    String tag = "";

    /**
     * Constructor. Place to initialise all histograms
     */
    public Histos() {
        this("");
    }

    /**
     * Constructor.
     * @param tag is used to ensure a unique identifier
     */
    public Histos(String tag) {

        this.tag = tag;

        h1_p = new H1F("reconstructed_p" + tag, "reconstructed p", 100, 0, 1000);
        h1_theta = new H1F("reconstructed_theta"  + tag, "reconstructed theta", 100, 0, 180);
        h1_phi = new H1F("reconstructed_phi"  + tag, "reconstructed phi", 100, 0, 360);
        h1_vz = new H1F("reconstructed_vz"  + tag, "reconstructed vz", 100, -24, 20);

        h1_p0 = new H1F("expected_p"  + tag, "expected p", 100, 0, 1000);
        h1_theta0 = new H1F("expected_theta"  + tag, "expected theta", 100, 0, 180);
        h1_phi0 = new H1F("expected_phi"  + tag, "expected phi", 100, 0, 360);
        h1_vz0 = new H1F("expected_vz"  + tag, "expected vz", 100, -35, 20);

        h1_delta_p = new H1F("delta_p"  + tag, "delta p", 100, -500, 500);
        h1_delta_theta = new H1F("delta_theta"  + tag, "delta theta", 100, -180, 180);
        h1_delta_phi = new H1F("delta_phi"  + tag, "delta phi", 100, -40, 40);
        h1_delta_vz = new H1F("delta_vz"  + tag, "delta vz", 100, -20, 10);

        h1_computing_time = new H1F("computing_time" + tag, "ALERTEngine computing time", 100, 0, 120);

        h1_residual = new H1F("residual" + tag, "residual", 100, -2.2, 2.2);

        /// --- title
        h1_computing_time.setTitleX("time (ms)");
        h1_computing_time.setTitleY("count");

        h1_residual.setTitleX("residual (mm)");
        h1_residual.setTitleY("count");

        h1_p.setTitleX("p (MeV)");
        h1_p.setTitleY("count");
        h1_theta.setTitleX("theta (deg)");
        h1_theta.setTitleY("count");
        h1_phi.setTitleX("phi (deg)");
        h1_phi.setTitleY("count");
        h1_vz.setTitleX("vz (cm)");
        h1_vz.setTitleY("count");

        h1_p0.setTitleX("p (MeV)");
        h1_p0.setTitleY("count");
        h1_theta0.setTitleX("theta (deg)");
        h1_theta0.setTitleY("count");
        h1_phi0.setTitleX("phi (deg)");
        h1_phi0.setTitleY("count");
        h1_vz0.setTitleX("vz (cm)");
        h1_vz0.setTitleY("count");

        h1_delta_p.setTitleX("delta p (MeV)");
        h1_delta_p.setTitleY("count");
        h1_delta_theta.setTitleX("delta theta (deg)");
        h1_delta_theta.setTitleY("count");
        h1_delta_phi.setTitleX("delta phi (deg)");
        h1_delta_phi.setTitleY("count");
        h1_delta_vz.setTitleX("delta vz (cm)");
        h1_delta_vz.setTitleY("count");

        /// --- Fill Map
        map_histo1d.put(h1_p.getName(), h1_p);
        map_histo1d.put(h1_theta.getName(), h1_theta);
        map_histo1d.put(h1_phi.getName(), h1_phi);
        map_histo1d.put(h1_vz.getName(), h1_vz);

        map_histo1d.put(h1_p0.getName(), h1_p0);
        map_histo1d.put(h1_theta0.getName(), h1_theta0);
        map_histo1d.put(h1_phi0.getName(), h1_phi0);
        map_histo1d.put(h1_vz0.getName(), h1_vz0);

        map_histo1d.put(h1_delta_p.getName(), h1_delta_p);
        map_histo1d.put(h1_delta_theta.getName(), h1_delta_theta);
        map_histo1d.put(h1_delta_phi.getName(), h1_delta_phi);
        map_histo1d.put(h1_delta_vz.getName(), h1_delta_vz);
        
        map_histo1d.put(h1_computing_time.getName(), h1_computing_time);
        map_histo1d.put(h1_residual.getName(), h1_residual);

    }

    /** 
     * Add the content of the incoming Histos object int the current one.
     * Very useful for the merge operation after the parallelisation.
     * @param histos
     */
    public void merge(Histos histos) {

        // this.h1_p.add(histos.h1_p);
        // this.h1_theta.add(histos.h1_theta);
        // this.h1_phi.add(histos.h1_phi);
        // this.h1_vz.add(histos.h1_vz);

        // this.h1_p0.add(histos.h1_p0);
        // this.h1_theta0.add(histos.h1_theta0);
        // this.h1_phi0.add(histos.h1_phi0);
        // this.h1_vz0.add(histos.h1_vz0);

        // this.h1_delta_p.add(histos.h1_delta_p);
        // this.h1_delta_theta.add(histos.h1_delta_theta);
        // this.h1_delta_phi.add(histos.h1_delta_phi);
        // this.h1_delta_vz.add(histos.h1_delta_vz);
        
        // this.h1_computing_time.add(histos.h1_computing_time);
        //this.h1_residual.add(histos.h1_residual);

        for (Map.Entry<String, H1F> entry : histos.map_histo1d.entrySet()) {
            String key = histos.removeTag(entry.getKey()) + tag;
            if (this.map_histo1d.containsKey(key)) {
                this.map_histo1d.get(key).add(entry.getValue());
            }
        }

    }

    // to be completed
    public void print() {
        System.out.println("Histos:");
        System.out.printf("Reconstructed  --->  p : %f MeV, theta : %f deg, phi : %f deg , vz : %f cm \n", h1_p.getMean(), h1_theta.getMean(), h1_phi.getMean(), h1_vz.getMean());
        System.out.printf("Expected       --->  p : %f MeV, theta : %f deg, phi : %f deg , vz : %f cm \n", h1_p0.getMean(), h1_theta0.getMean(), h1_phi0.getMean(), h1_vz0.getMean());
        System.out.printf("Delta          --->  Dp : %f MeV, Dtheta : %f deg, Dphi : %f deg , Dvz : %f cm \n", h1_delta_p.getMean(), h1_delta_theta.getMean(), h1_delta_phi.getMean(), h1_delta_vz.getMean());
        System.out.printf("event computing time : %f ms\n", h1_computing_time.getMean());
        System.out.printf("residuals : %f mm , std dev : %f mm\n", h1_residual.getMean(), h1_residual.getRMS());
        printListOfHistograms(true);
    }

    public void save(String outDir) throws IOException, DocumentException {

        // Renderer.width = 1500;
        // Renderer.height = 1200;

        for (Map.Entry<String, H1F> entry : map_histo1d.entrySet()) {
            String name = entry.getKey();
            H1F h = entry.getValue();
            Renderer.save_histogram(h, outDir + "/" + name, RendererOutputType.PNG);
        }

        // ArrayList<H1F> list = new ArrayList<>();
        // list.add(h1_p);
        // list.add(h1_p0);
        // Renderer.save_overlayed_histograms(list, "combined_histo_p", RendererOutputType.PNG);
        // Renderer.save_overlayed_histograms(list, "combined_histo_p", RendererOutputType.PDF);

    }

    public String getTag() {return tag;}

    public ArrayList<String> getListOfHistograms(boolean tag) {
        ArrayList<String> list = new ArrayList<>();
        for (String key : map_histo1d.keySet()) {
            if (tag) {
                list.add(removeTag(key));
            } else {
                list.add(key);
            }
        }
        return list;
    }

    /**
     * Print the list of histogram
     * @param tag ? do we remove tag for the names
     */
    public void printListOfHistograms(boolean tag) {
        System.out.printf("List of histograms (tag : %s)\n", this.tag);
        int counter = 0;
        for (String key : map_histo1d.keySet()) {
            counter++;
            if (tag) {
                System.out.printf("%2d)   %s\n", counter, removeTag(key));
            } else {
                System.out.printf("%2d)   %s\n", counter, key);
            }
        }
    }

    /**
     * Remove the {@link Histos#tag from a name}
     * @param name 
     * @return name without the tag (at the end), return the same String if it does not end with the tag
     */
    public String removeTag(String name) {
        if (name.endsWith(tag)) {
            return name.substring(0, name.length()-tag.length());
        } else {
            return name;
        }
    }

    /**
     * Access an histogram from the map by it key (no need to specify the tag)
     * @param name with or without the tag
     * @return correcponding H1F hitogram
     */
    public H1F getHistogram1D(String name) {
        return map_histo1d.get(removeTag(name) + tag);
    }
    
}
