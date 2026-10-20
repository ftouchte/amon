package io.github.ftouchte.kalmanFilter;

import java.util.ArrayList;
import java.util.HashMap;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.io.banks.AHDC_hits;
import org.jlab.io.banks.ATOF_hits;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.rec.ahdc.Hit.Hit;
import org.jlab.rec.ahdc.KalmanFilter.KFHit;
import org.jlab.rec.ahdc.KalmanFilter.KFitter;
import org.jlab.rec.ahdc.KalmanFilter.MaterialMap;
import org.jlab.rec.ahdc.KalmanFilter.Propagator;
import org.jlab.rec.ahdc.KalmanFilter.RadialKFHit;
import org.jlab.rec.ahdc.KalmanFilter.RungeKutta4;
import org.jlab.rec.alert.Track.Track;
import org.jlab.rec.alert.Track.TrackCandidate;
import org.jlab.service.ahdc.AHDCEngine;
import org.jlab.service.alert.ALERTEngine;

import java.awt.BasicStroke;
import java.awt.Color;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.util.ShapeUtils;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.clas.tracking.kalmanfilter.Material;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.axis.NumberAxis;


import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import java.io.File;
import java.io.IOException;

import javax.swing.*;
import java.awt.Shape;
import java.awt.geom.*;
import org.jfree.chart.util.ShapeUtils;


/**
 * This code is intended to be run on a filtered HIPO files
 * containing only elestic events
 * 
 * Just be fast...
 * 
 * @author Felix Touchte Codjo
 * @date October 13, 2026 
 * 
 * Monitor
 */

public class Monitor {


    public static void main(String[] args) throws Exception {
        String filename = "/work/clas12/users/touchte/data/kfmon/filtered_rec_clas_022712.evio.00040-00049.hipo";

        // Read hipo file
        HipoReader reader = new HipoReader();
        reader.open(filename);

        SchemaFactory factory = reader.getSchemaFactory();

        // Load the AHDC geometry
        DatabaseConstantProvider cp = new DatabaseConstantProvider(22712, "default");
        AlertDCDetector AHDC = (new AlertDCFactory()).createDetectorCLAS(cp);

        int nevents = 0;

        // Loop over events
        while (reader.hasNext()) {

            nevents++;
            if (nevents > 10) break;
                
            /// --- Load event and parse the factory
            Event event0 = new Event();
            reader.nextEvent(event0);
            DataEvent event = new HipoDataEvent(event0, factory);

            /// --- Reconstruct the Track and Hit objects : inspired by ALERTEngine
            DataBank trackBank = event.getBank("AHDC::track");
            DataBank hitBank = event.getBank("AHDC::hits");

            ArrayList<Track> AHDC_tracks = new ArrayList();

            for(int row = 0; row < trackBank.rows(); ++row) {
                int trackid = trackBank.getInt("trackid", row);
                ArrayList<Hit> AHDC_hits = new ArrayList();

                for(int hit_row = 0; hit_row < hitBank.rows(); ++hit_row) {
                    if (trackid == hitBank.getInt("trackid", hit_row)) {
                        int id = hitBank.getShort("id", hit_row);
                        int superlayer = hitBank.getByte("superlayer", hit_row);
                        int layer = hitBank.getByte("layer", hit_row);
                        int wire = hitBank.getInt("wire", hit_row);
                        int adc = hitBank.getInt("adc", hit_row);
                        double doca = hitBank.getDouble("doca", hit_row);
                        double time = hitBank.getDouble("time", hit_row);
                        double tot = hitBank.getDouble("timeOverThreshold", hit_row);
                        double gainCorr = 1.0;
                        Hit hit = new Hit(id, superlayer, layer, wire, doca, (double)adc / gainCorr, time);
                        hit.setWirePosition(AHDC);
                        hit.setTrackId(trackid);
                        hit.setADC((double)adc);
                        hit.setToT(tot);
                        AHDC_hits.add(hit);
                    }
                }

                TrackCandidate newCandidate = new TrackCandidate(AHDC_hits);
                Track newTrack = new Track(newCandidate);
                double[] vec = new double[]{(double)trackBank.getFloat("x", row), (double)trackBank.getFloat("y", row), (double)trackBank.getFloat("z", row), (double)trackBank.getFloat("px", row), (double)trackBank.getFloat("py", row), (double)trackBank.getFloat("pz", row)};
                newTrack.setPositionAndMomentumVec(vec);
                newTrack.set_trackId(trackid);
                AHDC_tracks.add(newTrack);
            } // end track and hit reconstrcution from hipo file



            /// --- Read the electron vz
            DataBank recBank = event.getBank("REC::Particle");
            double vertex = 0;
            double expected_p = 0;
            double expected_theta = 0;
            double expected_phi = 0;
            for (int row = 0; row < recBank.rows(); row++) {
                if (recBank.getInt("pid", row) == 11) {
                    vertex = 10.0 * recBank.getFloat("vz", row); // convert from cm to mm
                    double px = 1000*recBank.getFloat("px", row);
                    double py = 1000*recBank.getFloat("py", row);
                    double pz = 1000*recBank.getFloat("pz", row);
                    double p = Math.sqrt(px*px + py*py + pz*pz);
                    double theta = Math.acos(pz/p);
                    double phi = Math.atan2(py, px);

                    double beam_energy = 2239.49; // MeV
                    double px1 = -p*Math.sin(theta)*Math.cos(phi);
                    double py1 = -p*Math.sin(theta)*Math.sin(phi);
                    double pz1 = beam_energy - p*Math.cos(theta);

                    expected_p = Math.sqrt(px1*px1 + py1*py1 + pz1*pz1);
                    expected_theta = Math.acos(pz1/expected_p);
                    expected_phi = Math.atan2(py1, px1);

                    break;
                }
            }

            // Correct for the alignment
            vertex += 67;

            /// --- Run the Kalman Filter
            
            /// --- Some config objects
            DataBank runBank = event.getBank("RUN::config");
            HashMap<String, Material> materialHashMap = MaterialMap.generateMaterials();

            // Loop over tracks
            for (Track track : AHDC_tracks) {

                // initial fit value
                double x0 = (double)0.0F;
                double y0 = (double)0.0F;
                double z0 = vertex;
                double px0 = track.get_px();
                double py0 = track.get_py();
                double pz0 = track.get_pz();
                double[] y = new double[]{x0, y0, z0, px0, py0, pz0};


                // Setup the fitter
                PDGParticle proton = PDGDatabase.getParticleById(2212);
                double magfieldfactor = (double)runBank.getFloat("solenoid", 0);
                double magfield = (double)50.0F * magfieldfactor;
                double[] B = new double[]{(double)0.0F, (double)0.0F, magfield / (double)10.0F * 0.001};
                
                RungeKutta4 RK4 = new RungeKutta4(proton, 6, B);
                Propagator propagator = new Propagator(RK4);
                RealVector initialStateEstimate = new ArrayRealVector(y);
                // First 3 lines in mm^2; last 3 lines in MeV^2 (in the beamline)
                RealMatrix initialErrorCovariance = MatrixUtils.createRealMatrix(new double[][]{
                    {50  , 0.0 , 0.0 , 0.0 , 0.0 , 0.0}, 
                    {0.0 , 50  , 0.0 , 0.0 , 0.0 , 0.0}, 
                    {0.0 , 0.0 , 900 , 0.0 , 0.0 , 0.0}, 
                    {0.0 , 0.0 , 0.0 , 100 , 0.0 , 0.0}, 
                    {0.0 , 0.0 , 0.0 , 0.0 , 100 , 0.0}, 
                    {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 900}});
                track.setErrorCovarianceMatrix(initialErrorCovariance);
                KFitter TrackFitter = new KFitter(initialStateEstimate, initialErrorCovariance, propagator, materialHashMap);

                // Get AHDC hits
                ArrayList<Hit> AHDC_hits = track.getHits();

                // Generate beamline hit
                RadialKFHit hit_beam = new RadialKFHit(0, 0, vertex);
                double[] vz_error2 = {0.09, 1e10, 64}; // mm^2, radians^2, mm^2, error on r, phi, z
                RealMatrix measurementNoise = new Array2DRowRealMatrix(
                                                new double[][]{
                                                    {vz_error2[0], 0.0000      , 0.0000},
                                                    {0.0000      , vz_error2[1], 0.0000},
                                                    {0.0000      , 0.0000      , vz_error2[2]}
                                                });//3x3;
                hit_beam.setMeasurementNoise(measurementNoise);

                // Do KF iterations
                int Niter = 150;
                XYSeriesCollection prediction_dataset = new XYSeriesCollection();
                XYSeriesCollection correction_dataset = new XYSeriesCollection();
                XYSeriesCollection dataset = new XYSeriesCollection();
                
                XYSeriesCollection p_dataset = new XYSeriesCollection();
                XYSeriesCollection theta_dataset = new XYSeriesCollection();
                XYSeriesCollection phi_dataset = new XYSeriesCollection();

                XYSeries p_series = new XYSeries(String.format("p_series_evt_%d_trackid_%d", nevents, track.get_trackId()), false, true);
                XYSeries theta_series = new XYSeries(String.format("theta_series_evt_%d_trackid_%d", nevents, track.get_trackId()), false, true);
                XYSeries phi_series = new XYSeries(String.format("phi_series_evt_%d_trackid_%d", nevents, track.get_trackId()), false, true);

                for(int k = 0; k < Niter; ++k) {
                    
                    XYSeries predictions = new XYSeries(String.format("prediction_evt_%d_trackid_%d_iter_%d", nevents, track.get_trackId(), k), false, true);
                    XYSeries corrections = new XYSeries(String.format("correction_evt_%d_trackid_%d_iter_%d", nevents, track.get_trackId(), k), false, true);
                    XYSeries combined = new XYSeries(String.format("combined_evt_%d_trackid_%d_iter_%d", nevents, track.get_trackId(), k), false, true);

                    // Starting point
                    addZRPointToSeries(TrackFitter, combined);

                    //TrackFitter.ResetErrorCovariance(initialErrorCovariance.copy());



                    // Forward
                    for(Hit hit : AHDC_hits) {
                        TrackFitter.predict(hit, true);
                        addZRPointToSeries(TrackFitter, predictions);
                        addZRPointToSeries(TrackFitter, combined);
                        TrackFitter.correct(hit);
                        addZRPointToSeries(TrackFitter, corrections);
                        addZRPointToSeries(TrackFitter, combined);
                    }

                    // Backward
                    for(int i = AHDC_hits.size() - 2; i >= 0; i--) {
                        Hit hit = (Hit)AHDC_hits.get(i);
                        TrackFitter.predict(hit, false);
                        addZRPointToSeries(TrackFitter, predictions);
                        addZRPointToSeries(TrackFitter, combined);
                        TrackFitter.correct(hit);
                        addZRPointToSeries(TrackFitter, corrections);
                        addZRPointToSeries(TrackFitter, combined);
                    }

                    // Beamline
                    TrackFitter.predict(hit_beam, false);
                    addZRPointToSeries(TrackFitter, predictions);
                    addZRPointToSeries(TrackFitter, combined);
                    TrackFitter.correct(hit_beam);
                    addZRPointToSeries(TrackFitter, corrections);
                    addZRPointToSeries(TrackFitter, combined);

                    prediction_dataset.addSeries(predictions);
                    correction_dataset.addSeries(corrections);
                    dataset.addSeries(combined);

                    addMomemtumToSeries(TrackFitter, k, p_series, theta_series, phi_series);

                } // end loop over iter
                p_dataset.addSeries(p_series);
                theta_dataset.addSeries(theta_series);
                phi_dataset.addSeries(phi_series);

                // Do post fit propagation
                RealVector x_out = TrackFitter.getStateEstimationVector();
                track.setPositionAndMomentumVec(x_out.toArray());
                track.setErrorCovarianceMatrix(TrackFitter.getErrorCovarianceMatrix());
                KFitter PostFitPropagator = new KFitter(TrackFitter.getStateEstimationVector(), initialErrorCovariance, propagator, materialHashMap);

                for(Hit hit : AHDC_hits) {
                    PostFitPropagator.predict(hit, true);
                    hit.setResidual(PostFitPropagator.residual(hit));
                    hit.setResidual_LR(PostFitPropagator.residual_LR(hit));
                }

                // Expected track
                // to be done...

                String workDir = "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/kalmanFilter/output";
                String outName = String.format("%s/evt_%d_trackid_%d_ZR.png", workDir, nevents, track.get_trackId());
                create_zr_chart(prediction_dataset, correction_dataset, dataset, outName);

                // Momentum output
                String p_outName = String.format("%s/evt_%d_trackid_%d_p.png", workDir, nevents, track.get_trackId());
                create_niter_chart(p_dataset, "p (MeV)", expected_p, p_outName, Niter);

                String theta_outName = String.format("%s/evt_%d_trackid_%d_theta.png", workDir, nevents, track.get_trackId());
                create_niter_chart(theta_dataset, "theta (rad)", expected_theta, theta_outName, Niter);

                String phi_outName = String.format("%s/evt_%d_trackid_%d_phi.png", workDir, nevents, track.get_trackId());
                create_niter_chart(phi_dataset, "phi (rad)", expected_phi, phi_outName, Niter);

            } // end loop over tracks


        } // end loop over events


    } // end main

    public static void addZRPointToSeries (KFitter fitter, XYSeries series) {
        RealVector x_out = fitter.getStateEstimationVector();
        double x = x_out.getEntry(0);
        double y = x_out.getEntry(1);
        double z = x_out.getEntry(2);
        double r = Math.sqrt(x*x + y*y);
        series.add(z,r);
    }

    public static void addMomemtumToSeries (KFitter fitter, int Niter, XYSeries P, XYSeries Theta, XYSeries Phi) {
        RealVector x_out = fitter.getStateEstimationVector();
        double px = x_out.getEntry(3);
        double py = x_out.getEntry(4);
        double pz = x_out.getEntry(5);
        double p = Math.sqrt(px*px + py*py + pz*pz);
        double theta = Math.acos(pz/p);
        double phi = Math.atan2(py, px);
        P.add(Niter, p);
        Theta.add(Niter, theta);
        Phi.add(Niter, phi);
    }

    public static void create_zr_chart(XYSeriesCollection prediction_dataset, XYSeriesCollection correction_dataset, XYSeriesCollection dataset, String outName) {
        // Create the chart with the first dataset
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Track fit",       // title
                "z",               // x-axis label
                "r",               // y-axis label
                dataset,           // dataset 0 (combined)
                PlotOrientation.VERTICAL,
                true, true, false);

        XYPlot plot = chart.getXYPlot();

        // The last serie is the last to be plotted
        plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);

        //Axis limit
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();

        xAxis.setRange(-200.0, 200.0);   // xlim
        yAxis.setRange(0.0, 80.0);   // ylim




        // Attach the other datasets
        plot.setDataset(1, prediction_dataset);
        plot.setDataset(2, correction_dataset);

        XYLineAndShapeRenderer r0 = new XYLineAndShapeRenderer(true, false);
        XYLineAndShapeRenderer r1 = new XYLineAndShapeRenderer(false, true);
        XYLineAndShapeRenderer r2 = new XYLineAndShapeRenderer(false, true);

        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            int n = dataset.getSeriesCount();
            float t = (n > 1) ? i / (float) (n - 1) : 1f;
            Color c = Color.getHSBColor(0.6f, 0.2f + 0.8f * t, 0.5f + 0.5f * t);
            r0.setSeriesPaint(i, c);   // (indice, couleur)
            r0.setSeriesStroke(i, new BasicStroke(1.0f));
        }
        for (int i = 0; i < prediction_dataset.getSeriesCount(); i++) {
            r1.setSeriesPaint(i, java.awt.Color.BLUE);
            double s = 3.0;
            r1.setSeriesShape(i, ShapeUtils.createUpTriangle((float) s));
        }
        for (int i = 0; i < correction_dataset.getSeriesCount(); i++) {
            r2.setSeriesPaint(i, java.awt.Color.RED);
            double s = 2.0;
            r2.setSeriesShape(i, new Rectangle2D.Double(-s, -s, 2*s, 2*s));
        }

        // remove legend
        for (int i = 1; i < dataset.getSeriesCount() - 1; i++) {
            r0.setSeriesVisibleInLegend(i, false);
            r1.setSeriesVisibleInLegend(i, false);
            r2.setSeriesVisibleInLegend(i, false);
        }

        plot.setRenderer(0, r0);
        plot.setRenderer(1, r1);
        plot.setRenderer(2, r2);

        try {
            ChartUtils.saveChartAsPNG(new File(outName), chart, 1200, 800);
        } catch (IOException e) {
            e.printStackTrace();
        }

    } // end create zr chart

    public static void create_niter_chart(XYSeriesCollection dataset, String valName, double expected_value, String outName, int max_Niter) {
        // Create the chart with the first dataset
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Track fit",       // title
                "Niter",               // x-axis label
                valName,               // y-axis label
                dataset,           // dataset 0 (combined)
                PlotOrientation.VERTICAL,
                true, true, false);

        XYPlot plot = chart.getXYPlot();

        // The last serie is the last to be plotted
        plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);

        //Axis limit
        // NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        // NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();

        // xAxis.setRange(-200.0, 200.0);   // xlim
        // yAxis.setRange(0.0, 80.0);   // ylim

        XYSeriesCollection line = new XYSeriesCollection();
        XYSeries line_Series = new XYSeries("expected");
        line_Series.add(0, expected_value);
        line_Series.add(max_Niter, expected_value);
        line.addSeries(line_Series);

        // Attach the other datasets
        plot.setDataset(1, line);

        XYLineAndShapeRenderer r0 = new XYLineAndShapeRenderer(true, true);
        XYLineAndShapeRenderer r1 = new XYLineAndShapeRenderer(true, false);

        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            int n = dataset.getSeriesCount();
            float t = (n > 1) ? i / (float) (n - 1) : 1f;
            Color c = Color.getHSBColor(0.6f, 0.2f + 0.8f * t, 0.5f + 0.5f * t);
            r0.setSeriesPaint(i, c);   // (indice, couleur)
            r0.setSeriesStroke(i, new BasicStroke(1.0f));
            double s = 2.0;
            r0.setSeriesShape(i, new Rectangle2D.Double(-s, -s, 2*s, 2*s));
        }
        for (int i = 0; i < line.getSeriesCount(); i++) {
            r1.setSeriesPaint(i, java.awt.Color.RED);
            r1.setSeriesStroke(i, new BasicStroke(1.0f));
        }

        // remove legend
        for (int i = 1; i < dataset.getSeriesCount() - 1; i++) {
            r0.setSeriesVisibleInLegend(i, false);
            r1.setSeriesVisibleInLegend(i, false);
        }

        plot.setRenderer(0, r0);
        plot.setRenderer(1, r1);
        try {
            ChartUtils.saveChartAsPNG(new File(outName), chart, 1200, 800);
        } catch (IOException e) {
            e.printStackTrace();
        }

    } // end create chart
}
