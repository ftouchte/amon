package io.github.ftouchte.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Supplier;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.service.ahdc.AHDCEngine;
import org.jlab.service.alert.ALERTEngine;

import io.github.ftouchte.alignment.ResultsOverIterations.CCDB_TYPE;
import io.github.ftouchte.filtering.AlertElasticAnalyser;
import io.github.ftouchte.filtering.AlertGoodTrackFinder;
import io.github.ftouchte.filtering.AlertTrackSelector;
import io.github.ftouchte.fitting.CrystalBall;
import io.github.ftouchte.fitting.CrystallBallFitter;
import io.github.ftouchte.utils.ParticleRow;
import io.github.ftouchte.utils.Units;
import io.github.ftouchte.utils.fOptions;


/**
 * This code is used in parallel of coatjava{branch: ahdc/alignment}. We use coatjva to access the 
 * AHDC geometry and the revisited processEvent() method of some classes.
 * 
 * @note This class contains analysis code ready to used for:
 * - layer alignment
 * - wire alignment
 * - position scan
 */
public class AhdcAlignmentAnalyser {

    /**
     * This the main method. It is used to run an iteration of the alignment procedure.
     * This code uses parallelism to reduce the computing time. The user is invited to
     * manage the parameters : {@link #nThreads} and {@link #queue_capacity}
     * 
     * @param niter Iteration number
     * @param results {@link ResultsOverIterations} where to store data over iteration
     * @param inFiles HIPO files containing all events to be processed
     * @param outDir where to store plots
     * @param clas_alignment position of the center of CLAS with respect to the center of ALERT. Useful for {@link #scan_ahdc_position(String[])}
     * @param flag_do_fit perform a KF fit or do a simple propagation
     * @param trackSelectorFactory {@link AlertTrackSelector}
     * @param analysers {@link AhdcAlignmentAnalyser.HistoAnalyser}
     */
    static void run(int niter, ResultsOverIterations results, ArrayList<String> inFiles, String outDir, double clas_alignment, boolean flag_do_fit, Supplier<AlertTrackSelector> trackSelectorFactory, HistoAnalyser... analysers) {
        
        /// --- Parallelizer
        BlockingQueue<DataEvent> queue = new ArrayBlockingQueue<>(queue_capacity);
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);

        DataEvent EVT_POISON = new HipoDataEvent(new Event());

        /// --- Worker
        // Create as many future (thread result) as the number of threads
        List<Future<Histos>> futures = new ArrayList<>();

        for (int i = 0; i < nThreads; i++) {
            futures.add(pool.submit(() -> {
                Histos local_histos = new Histos(niter);
                try {    
                    System.out.println("\033[1;32m **** Create new thread : " + Thread.currentThread().getName() + " \033[0m");
                    AlertTrackSelector trackSelector = trackSelectorFactory.get();

                    // Engines
                    AHDCEngine ahdcEngine = new AHDCEngine();
                    ALERTEngine alertEngine = new ALERTEngine();

                    // init ccdb and geometry
                    synchronized (CCDB_LOCK) {
                        ahdcEngine.init();
                        ahdcEngine.detectorChanged(22712);
                        alertEngine.init();
                        alertEngine.detectorChanged(22712);
                    }

                    // ALERT engine
                    alertEngine.set_clas_alignment(clas_alignment);
                    alertEngine.setStepSize(stepSize);
                    alertEngine.set_KF_Niter(25);

                    int nevents = 0;

                    while (true) {
                        DataEvent event = queue.take();

                        nevents++;
                        if (nevents % frequency_of_event_login == 0) {
                            System.out.println("\033[1;32m > " + Thread.currentThread().getName() + " : \033[0m" + nevents + " events");
                        }
                    
                        if (event == EVT_POISON) break;

                        fill_histos(local_histos, event, trackSelector, alertEngine, ahdcEngine, flag_do_fit);

                    }
                } catch (Throwable t) {
                    System.err.println("\033[1;31m WORKER EXCEPTION/ERROR: \033[0m");
                    t.printStackTrace();
                    throw t; // Throwable peut être re-lancé directement
                }

                return local_histos;
            }));
        }

        int nevents = 0;

        /// --- Producer
        for (String file : inFiles) {
            // Initialise HIPO reader for the file
            HipoReader reader = new HipoReader();
            reader.open(file);
            /// Read bank definitions
            SchemaFactory factory = reader.getSchemaFactory();


            while (reader.hasNext() && nevents < nevents_max) {
                // load event
                nevents++;
                Event event0 = new Event();
                reader.nextEvent(event0);
                DataEvent event = new HipoDataEvent(event0, factory);

                // put event in the queue
                try {
                    queue.put(event);
                } catch (InterruptedException e) {
                    System.out.println("\033[1;31m Error: fail to put event in the queue \033[0m");
                    e.printStackTrace();
                }
            } // end loop over over for this files
        } // end loop over files and events

        /// --- Stop worker with EVT_POISON signals
        for (int i = 0; i < nThreads; i++) {
            try {
                queue.put(EVT_POISON);
            } catch (InterruptedException e) {
                System.out.println("\033[1;31m Error: fail to put EVT_POISON in the queue \033[0m");
                e.printStackTrace();
            }
        }

        /// --- Merge histograms
        Histos global_histos = new Histos(niter);
        for (Future<Histos> f : futures) {
            try {
                Histos local_histos = f.get();
                global_histos.merge(local_histos);
            } catch (InterruptedException e) {
                System.out.println("\033[1;31m Error: fail to put EVT_POISON in the queue \033[0m");
                e.printStackTrace();
            } catch (ExecutionException e) {
                System.out.println("\033[1;31m Worker failed! \033[0m");
                e.getCause().printStackTrace(); 
            }   
        }

        /// --- Stop parallelizer
        pool.shutdown();
        try {
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS); // attend que tous les threads finissent
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // /// --- Analyse global histograms
        // analyse_global_histograms(global_histos, outDir, niter, results);

        // /// --- Analyse histograms per layer
        // analyse_layer_histograms(global_histos, outDir, niter, results);
        // analyse_layer_2D_histograms(global_histos, outDir, niter, results);

        // /// --- Analyse histograms per wire
        // if (flag_run_wire_alignment) {
        //     analyse_wire_histograms(global_histos, outDir, niter, results);
        //     //analyse_wire_2D_histograms(global_histos, outDir, niter, results);
        // }

        /// --- Run whichever analyses were requested
        for (HistoAnalyser analyser : analysers) {
            analyser.analyse(global_histos, outDir, niter, results);
        }
        
       
    }

    static void analyse_global_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        check_output_dir(outDir + "/global/");
        String globalOutDir = outDir + "/global/iter-" + niter;
        check_output_dir(globalOutDir);

        /// track theta
        {
            H1F h = global_histos.h1_track_theta;
            double mean = h.getMean();
            double width = h.getRMS();
            double clas_alignment = results.clas_alignment;
            h.setTitle(String.format("mean : %.5f, width : %.5f, clas alignment : %.1f", mean, width, -0.1*clas_alignment));
            save_histo1D(h, null,  globalOutDir + "/track_theta.pdf");
        }

        /// track delta phi
        {
            H1F h = global_histos.h1_track_delta_phi;
            double mean = h.getMean();
            double width = h.getRMS();
            h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            save_histo1D(h, null,  globalOutDir + "/track_delta_phi.pdf");
        }

        /// delta vz
        {
            H1F h = global_histos.h1_delta_vz;
            double mean = h.getMean();
            double width = h.getRMS();
            h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            save_histo1D(h, null,  globalOutDir + "/track_delta_vz.pdf");

            results.mean_delta_vz = mean;
            results.width_delta_vz = width;
        }

    }

    static void analyse_layer_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        check_output_dir(outDir + "/layers/");
        String layerOutDir = outDir + "/layers/iter-" + niter;
        check_output_dir(layerOutDir);

        /// --- residual LR
        EmbeddedCanvas c = new EmbeddedCanvas(1500, 1200);
        c.divide(3, 3);
        GraphErrors g_layer = new GraphErrors();
        GraphErrors g_layer_cost = new GraphErrors();
        for (int i = 0; i < global_histos.h1_residual_LR_per_layers.size(); i++) {

            H1F h = global_histos.h1_residual_LR_per_layers.get(i);
            
            double mean = 0;
            double width = -1;
            GraphErrors gr = null;
            Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "residual LR, layer " + i + "/8 in analyse_layer_histograms()...");
            if (fitResult != null) {
                CrystalBall cb = fitResult.getFirst();
                // System.out.println("* Fit result : residual LR layer " + i);
                // cb.print();
                gr = fitResult.getSecond();
                mean = cb.getMu();
                width = cb.getSigma();

                g_layer.addPoint(i, mean, 0, 0);
                g_layer_cost.addPoint(i, cb.getFitCost(), 0, 0);
            }

            // store residual for the current rotation angle 
            if (i > 0) {                      
                results.layer_residuals[i-1] = mean;
                h.setTitle(String.format("mean : %.5f, width : %.5f, alpha : %.4f -> %.4f deg", mean, width, results.layer_angles_start[i-1], results.layer_angles_end[i-1]));
            }
            if (i == 0) {
                h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            }           

            // Draw histograms
            c.cd(i);
            c.draw(h);
            if (gr != null)
                c.draw(gr, "same L");
            
        }
        c.save(layerOutDir + "/residual_LR_iter_" + niter + ".pdf");
        System.out.println("residual_LR_iter_" + niter + ".pdf created");

        EmbeddedCanvas c_layer = new EmbeddedCanvas(1200, 900);
        g_layer.setTitleX("layer");
        g_layer.setTitleY("mean residual LR");
        g_layer.setTitle("Mean residual LR versus layer");
        c_layer.draw(g_layer);
        c_layer.save(layerOutDir + "/monitoring_residual_LR_means_iter_" + niter + ".pdf");
        System.out.println(layerOutDir + "/monitoring_residual_LR_mean_bilan_iter_" + niter + ".pdf");

        EmbeddedCanvas c_layer_cost = new EmbeddedCanvas(1200, 900);
        g_layer_cost.setTitleX("layer");
        g_layer_cost.setTitleY("cost");
        g_layer_cost.setTitle("Least squares fit cost versus layer");
        c_layer_cost.draw(g_layer_cost);
        c_layer_cost.save(layerOutDir + "/monitoring_fit_quality_residual_LR_iter_" + niter + ".pdf");
        System.out.println(layerOutDir + "/monitoring_fit_quality_residual_LR_iter_" + niter + ".pdf");

        /// residual normal
        EmbeddedCanvas c2 = new EmbeddedCanvas(1500, 1200);
        c2.divide(3, 3);
        for (int i = 0; i < global_histos.h1_residual_per_layers.size(); i++) {
            H1F h = global_histos.h1_residual_per_layers.get(i);
            // fit
            double mean = h.getMean();
            double width = h.getRMS();
            GraphErrors gr = null;
            Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "residual normal, layer " + i + "/8 in analyse_layer_histograms()...");
            if (fitResult != null) {
                CrystalBall cb = fitResult.getFirst();
                // System.out.println("* Fit result : residual normal layer " + i);
                // cb.print();
                gr = fitResult.getSecond();
                mean = cb.getMu();
                width = cb.getSigma();
            }
            // title
            h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));
            // draw
            c2.cd(i);
            c2.draw(h);
            c2.draw(gr, "same L");
        }
        c2.save(layerOutDir + "/residual_normal_iter_" + niter + ".pdf");
        System.out.println("residual_normal_iter_" + niter + ".pdf created");

    }

    static void analyse_wire_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {

        // create a new dir
        check_output_dir(outDir + "/wires/");
        String wireOutDir = outDir + "/wires/iter-" + niter;
        check_output_dir(wireOutDir);

        /// --- residual LR
        // (4x3)x48 = 576
        ArrayList<EmbeddedCanvas> canvas = new ArrayList<>(); // 48 canvas
        for (int i = 0; i < 48; i++) {
            EmbeddedCanvas c0 = new EmbeddedCanvas(1500, 1600);
            c0.divide(3, 4);
            canvas.add(c0);
        }
        
        GraphErrors g_wire = new GraphErrors();
        GraphErrors g_wire_cost = new GraphErrors();

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(wireOutDir + "/bilan_residual_LR_iter_" + niter + ".csv"))) {
            writer.write("# residual LR per wires");
            writer.write("# tilte (num, sector, layer, component, mean, width, integral, cost");
            writer.newLine();
            writer.newLine();
            for (int i = 0; i < global_histos.h1_residual_LR_per_wires.size(); i++) {
                H1F h = global_histos.h1_residual_LR_per_wires.get(i);

                AhdcWireId identifier = new AhdcWireId(i);
                
                String line = String.format("%3d, %2d, %2d, %2d,   ", i, identifier.sector, identifier.layer, identifier.component);
                
                double mean = 0; // do nothing if the fit failed
                double width = -1;
                double cost = -1;
                GraphErrors gr = null;
                if (h.integral() > 1000) {
                    Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "residual LR, wire " + (i+1) + "/576 in analyse_wire_histograms()...");
                    if (fitResult != null) {
                        CrystalBall cb = fitResult.getFirst();
                        // System.out.println("* Fit result : residual LR layer " + i);
                        // cb.print();
                        gr = fitResult.getSecond();
                        mean = cb.getMu();
                        width = cb.getSigma();
                        cost = cb.getFitCost();

                        g_wire.addPoint(i, mean, 0, 0);
                        g_wire_cost.addPoint(i, cb.getFitCost(), 0, 0);
                    }
                }
                // store residual for the current rotation angle                  
                results.wire_residuals[i] = mean;
                h.setTitle(String.format("mean : %.5f, width : %.5f, Dalpha : %.4f deg", mean, width, results.wire_angles[i]));
                
                line += String.format("%f, %f, %f, %f", mean, width, h.integral(), cost);
                writer.write(line);
                writer.newLine();

                // draw historgrams
                int canvas_num = i / 12; // between 0 and 47
                int canvas_frame = i % 12; // bbetween 0 and 11
                canvas.get(canvas_num).cd(canvas_frame);
                canvas.get(canvas_num).draw(h);
                if (gr != null)
                    canvas.get(canvas_num).draw(gr, "same L");
            }
        } catch (IOException e) {
            System.err.println("Failure : save residual LR fit results");
            e.printStackTrace();
        }

        // save histograms
        for (int i = 0; i < 48; i++) {
            canvas.get(i).save(wireOutDir + String.format("/resilual_LR_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
            System.out.println(String.format("resilual_LR_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
        }

        EmbeddedCanvas c_wire = new EmbeddedCanvas(1200, 900);
        g_wire.setTitleX("wire");
        g_wire.setTitleY("mean residual LR");
        g_wire.setTitle("Mean residual LR versus wire");
        c_wire.draw(g_wire);
        c_wire.save(wireOutDir + "/summary-residual-LR-wire-iter-" + niter + ".pdf");
        System.out.println(wireOutDir + "/summary-residual-LR-wire-iter-" + niter + ".pdf");

        EmbeddedCanvas c_wire_cost = new EmbeddedCanvas(1200, 900);
        g_wire_cost.setTitleX("wire");
        g_wire_cost.setTitleY("cost");
        g_wire_cost.setTitle("Least squares fit cost versus wire");
        c_wire_cost.draw(g_wire_cost);
        c_wire_cost.save(wireOutDir + "/monitoring_fit_quality_residual_LR_iter_" + niter + ".pdf");
        System.out.println(wireOutDir + "/monitoring_fit_quality_residual_LR_iter_" + niter + ".pdf");
    }

    static void analyse_layer_2D_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        // same outdir as for 1D histograms
        check_output_dir(outDir + "/layers/");
        String layerOutDir = outDir + "/layers/iter-" + niter;
        check_output_dir(layerOutDir);
        String slicesDir = layerOutDir + "/corr-slices";
        check_output_dir(slicesDir);

        /// --- correlation vz residual LR
        EmbeddedCanvas c = new EmbeddedCanvas(1500, 1200);
        c.divide(3, 3);

        for (int i = 0; i < global_histos.h2_corr_vz_residual_LR_per_layers.size(); i++) {
            H2F h2_initial = global_histos.h2_corr_vz_residual_LR_per_layers.get(i);
            H2F h2 = h2_initial.rebinX(5); // regroup x axis by group of 5 bins
            
            GraphErrors gr = new GraphErrors();
            gr.setMarkerSize(4);
            gr.setLineColor(1); // black

            GraphErrors gr_bis = new GraphErrors();
            gr_bis.setMarkerSize(4);
            gr_bis.setLineColor(1); // black


            for (int bin = 0; bin < h2.getXAxis().getNBins(); bin++) {
                
                H1F h = h2.sliceX(bin);

                // fit
                double mean = h.getMean();
                double width = h.getRMS();
                GraphErrors func = null;
                Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "correlation vz and residual LR, layer " + i + "/8, slice " + bin + " in analyse_wire_histograms()...");
                if (fitResult != null) {
                    CrystalBall cb = fitResult.getFirst();
                    // System.out.println("* Fit result : corr vz residual LR layer + " + i + ",  bin " + bin);
                    // cb.print();
                    func = fitResult.getSecond();
                    mean = cb.getMu();
                    width = cb.getSigma();
                }
                h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));

                // store results
                //if (width < 0.5){
                if (true){
                    gr.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, width);
                    gr_bis.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, 0);
                }
                h.setTitleX("residual LR (mm)");
                h.setTitleY("count");

                // output
                save_histo1D(h, func, String.format("%s/residual-LR-layer-%d-slice-%d.pdf", slicesDir, i > 0 ? AhdcWireId.number2layer(i) : 0, bin));
            } // end loop over bin

            // now fit the graph with a strait line
            WeightedObservedPoints observedPoints = new WeightedObservedPoints();
            for (int pt = 0; pt < gr_bis.getDataSize(0); pt++) {
                double x = gr_bis.getDataX(pt);
                double y = gr_bis.getDataY(pt);
                if (x > -130 && x < 100)
                    observedPoints.add(x, y);
            }
            PolynomialCurveFitter polFitter = PolynomialCurveFitter.create(1);
            double[] params = polFitter.fit(observedPoints.toList());

            F1D func = new F1D("strait-fit", String.format("%f + %f*x", params[0], params[1]), -215, 160);
            func.setLineColor(2);
            func.setLineWidth(4);
            
            // draw
            c.cd(i);
            c.getPad(i).setPalette("kBird");
            c.getPad(i).getAxisFrame().setDrawAxisZ(false);
            h2_initial.setTitle(String.format("300*slope : %f, constant : %f", params[1]*300, params[0]));
            c.draw(h2_initial);
            c.draw(func, "same L");
            c.draw(gr_bis, "same P");
            c.getPad(i).getAxisY().setRange(-0.5, 0.5);
            c.getPad(i).getAxisFrame().setDrawAxisZ(false); // again

            // residual analysis based on the fit
            if (i > 0) {
                double slope = params[1];
                double constant = params[0];
                double z_min = -188;
                double z_max = 162.5;
                double residual_start = slope*z_min + constant; 
                double residual_end = slope*z_max + constant; 

                results.layer_residuals_start[i-1] = residual_start;
                results.layer_residuals_end[i-1] = residual_end;

                results.layer_residuals_slope[i-1] = slope;
                results.layer_residuals_constant[i-1] = constant;
            }

        } // end loop over h2

        c.save(layerOutDir + "/corr_vz_residual_LR_iter_" + niter + ".pdf");
        System.out.println("corr_vz_residual_LR_iter_" + niter + ".pdf created");

        /// --- time2distance
        EmbeddedCanvas c2 = new EmbeddedCanvas(1500, 1200);
        c2.divide(3, 3);
        ArrayList<double[]> t2d_params_list = new ArrayList<>();
        for (int i = 0; i < global_histos.h2_time2distance_per_layers.size(); i++) {
            H2F h2_initial = global_histos.h2_time2distance_per_layers.get(i);
            H2F h2 = h2_initial.rebinX(5); // regroup x axis by group of 5 bins

            GraphErrors gr = new GraphErrors();
            gr.setMarkerSize(4);
            gr.setLineColor(1); // black

            GraphErrors gr_bis = new GraphErrors();
            gr_bis.setMarkerSize(4);
            gr_bis.setLineColor(1); // black

            for (int bin = 0; bin < h2.getXAxis().getNBins(); bin++) {
                
                H1F h = h2.sliceX(bin);
                h.setTitleX("distance (mm)");
                h.setTitleY("count");

                double mean = h.getMean();
                double width = h.getRMS();
                GraphErrors func = null;
                Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "time2distance, layer " + i + "/8, slice " + bin + " in analyse_wire_histograms()...");
                if (fitResult != null) {
                    CrystalBall cb = fitResult.getFirst();
                    // System.out.println("* Fit result : time2distance layer + " + i + ",  bin " + bin);
                    // cb.print();
                    func = fitResult.getSecond();
                    mean = cb.getMu();
                    width = cb.getSigma();
                }

                // store results
                //if (width < 0.5){
                if (true){
                    gr.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, width);
                    gr_bis.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, 0);
                }
                h.setTitle(String.format("mean : %f , width %f", mean, width));

                // output
                save_histo1D(h, func, String.format("%s/time2distance-layer-%d-slice-%d.pdf", slicesDir, i > 0 ? AhdcWireId.number2layer(i) : 0, bin));
            } // end loop over bin

            // Data point to be fit
            ArrayList<Double> xList = new ArrayList<>();
            ArrayList<Double> yList = new ArrayList<>();
            // xList.add(0.0);
            // yList.add(0.0);
            for (int pt = 0; pt < gr_bis.getDataSize(0)-2; pt++) { // exclude the 2 end points
                double x = gr_bis.getDataX(pt);
                double y = gr_bis.getDataY(pt);
                if (x < 230) { // filtering
                    xList.add(x);
                    yList.add(y);
                }
            }
            double[] xData = new double[xList.size()];
            double[] yData = new double[xList.size()];

            for (int bin = 0; bin < xData.length; bin++) {
                xData[bin] = xList.get(bin);
                yData[bin] = yList.get(bin);
            }

            // Fit
            GraphErrors func = null;
            double[] initialValues = results.time2distance.get(AhdcWireId.number2layer(Math.max(i, 1))); // use the t2d of the wire 1 in the corresponding layer as initialization; if i = 0 (integrated over all wire), initialize with layer 11

            try {
                double[] fittedParameters = fit_time2distance(xData, yData, initialValues.clone());
                int Npts = 1000;
                double tmin = 0;
                double tmax = 250;
                func = new GraphErrors();
                func.setLineColor(2);
                func.setLineThickness(1);
                func.setMarkerSize(0);
                for (int bin = 0; bin < Npts; bin++) {
                    double time = tmin + (tmax-tmin)*bin/(Npts-1);
                    func.addPoint(time, eval_t2d(time, fittedParameters), 0, 0);
                }
                t2d_params_list.add(fittedParameters);
            } catch (Exception e) {
                System.out.println("Failed to fit time2disatnce: layer " + i);
                // System.out.println("Error: " + e.getMessage());
                // e.printStackTrace();
                t2d_params_list.add(initialValues.clone());
            }
            
            c2.cd(i);
            c2.getPad(i).setPalette("kBird");
            c2.getPad(i).getAxisFrame().setDrawAxisZ(false);
            h2_initial.setTitle("time2distance");
            c2.draw(h2_initial);
            if (func != null) c2.draw(func, "same L");
            c2.draw(gr_bis, "same P");
            c2.getPad(i).getAxisY().setRange(0, 3);
            c2.getPad(i).getAxisX().setRange(0, 250);
            c2.getPad(i).getAxisFrame().setDrawAxisZ(false); // again

        }
        c2.save(layerOutDir + "/time2distance_iter_" + niter + ".pdf");
        System.out.println("time2distance_LR_iter_" + niter + ".pdf created");
        save_time2distance(t2d_params_list, results);

    }

    static void analyse_wire_2D_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        // same outdir as for 1D histograms
        check_output_dir(outDir + "/wires/");
        String wireOutDir = outDir + "/wires/iter-" + niter;
        check_output_dir(wireOutDir);
        String slicesDir = wireOutDir + "/corr-slices";
        check_output_dir(slicesDir);

        /// --- correlation vz residual LR
        // (4x3)x48 = 576
        ArrayList<EmbeddedCanvas> canvas = new ArrayList<>(); // 48 canvas
        for (int i = 0; i < 48; i++) {
            EmbeddedCanvas c0 = new EmbeddedCanvas(1500, 1600);
            c0.divide(3, 4);
            canvas.add(c0);
        }

        for (int i = 0; i < global_histos.h2_corr_vz_residual_LR_per_wires.size(); i++) {
            H2F h2_initial = global_histos.h2_corr_vz_residual_LR_per_wires.get(i);
            H2F h2 = h2_initial.rebinX(5); // regroup x axis by group of 5 bins
            
            GraphErrors gr = new GraphErrors();
            gr.setMarkerSize(4);
            gr.setLineColor(1); // black

            GraphErrors gr_bis = new GraphErrors();
            gr_bis.setMarkerSize(4);
            gr_bis.setLineColor(1); // black

            AhdcWireId identifier = new AhdcWireId(i);

            for (int bin = 0; bin < h2.getXAxis().getNBins(); bin++) {
                
                H1F h = h2.sliceX(bin);

                // fit
                double mean = h.getMean();
                double width = h.getRMS();
                GraphErrors func = null;
                Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "correlation vz and residual LR, wire " + (i+1) + "/576, slice " + bin + " in analyse_wire_histograms()...");
                if (fitResult != null) {
                    CrystalBall cb = fitResult.getFirst();
                    // System.out.println("* Fit result : corr vz residual LR layer + " + i + ",  bin " + bin);
                    // cb.print();
                    func = fitResult.getSecond();
                    mean = cb.getMu();
                    width = cb.getSigma();
                }
                h.setTitle(String.format("mean : %.5f, width : %.5f", mean, width));

                // store results
                //if (width < 0.5){
                if (true){
                    gr.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, width);
                    gr_bis.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, 0);
                }
                h.setTitleX("residual LR (mm)");
                h.setTitleY("count");

                // output
                save_histo1D(h, func, String.format("%s/residual-LR-wire-L%dW%d-slice-%d.pdf", slicesDir, identifier.layer, identifier.component, bin));
            } // end loop over bin

            // now fit the graph with a strait line
            WeightedObservedPoints observedPoints = new WeightedObservedPoints();
            for (int pt = 2; pt < gr_bis.getDataSize(0)-2; pt++) { // the 4 end points are excluded
                observedPoints.add(gr_bis.getDataX(pt), gr_bis.getDataY(pt));
            }
            PolynomialCurveFitter polFitter = PolynomialCurveFitter.create(1);
            double[] params = polFitter.fit(observedPoints.toList());

            F1D func = new F1D("strait-fit", String.format("%f + %f*x", params[0], params[1]), -215, 160);
            func.setLineColor(2);
            func.setLineWidth(4);
            
            // draw
            int canvas_num = i / 12; // between 0 and 47
            int canvas_frame = i % 12; // between 0 and 11
            canvas.get(canvas_num).cd(canvas_frame);
            canvas.get(canvas_num).getPad(canvas_frame).setPalette("kBird");
            canvas.get(canvas_num).getPad(canvas_frame).getAxisFrame().setDrawAxisZ(false);
            h2_initial.setTitle(String.format("300*slope : %f, constant : %f", params[1]*300, params[0]));
            canvas.get(canvas_num).draw(h2_initial);
            canvas.get(canvas_num).draw(func, "same L");
            canvas.get(canvas_num).draw(gr_bis, "same P");
            canvas.get(canvas_num).getPad(canvas_frame).getAxisY().setRange(-0.5, 0.5);
            canvas.get(canvas_num).getPad(canvas_frame).getAxisFrame().setDrawAxisZ(false); // again

            // residual analysis based on the fit
            double slope = params[1];
            double constant = params[0];
            double z_min = -188;
            double z_max = 162.5;
            double residual_start = slope*z_min + constant; 
            double residual_end = slope*z_max + constant; 

            results.wire_residuals_start[i] = residual_start;
            results.wire_residuals_end[i] = residual_end;

            results.wire_residuals_slope[i] = slope;
            results.wire_residuals_constant[i] = constant;

        } // end loop over h2

        // save histograms
        for (int i = 0; i < 48; i++) {
            canvas.get(i).save(wireOutDir + String.format("/corr_vz_resilual_LR_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
            System.out.println(String.format("corr_vz_resilual_LR_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
        }

        /// --- time2distance
        // (4x3)x48 = 576
        ArrayList<EmbeddedCanvas> canvas2 = new ArrayList<>(); // 48 canvas
        for (int i = 0; i < 48; i++) {
            EmbeddedCanvas c0 = new EmbeddedCanvas(1500, 1600);
            c0.divide(3, 4);
            canvas2.add(c0);
        }
        ArrayList<double[]> t2d_params_list = new ArrayList<>();
        for (int i = 0; i < global_histos.h2_time2distance_per_wires.size(); i++) {
            H2F h2_initial = global_histos.h2_time2distance_per_wires.get(i);
            H2F h2 = h2_initial.rebinX(5); // regroup x axis by group of 5 bins

            GraphErrors gr = new GraphErrors();
            gr.setMarkerSize(4);
            gr.setLineColor(1); // black

            GraphErrors gr_bis = new GraphErrors();
            gr_bis.setMarkerSize(4);
            gr_bis.setLineColor(1); // black

            AhdcWireId identifier = new AhdcWireId(i);

            for (int bin = 0; bin < h2.getXAxis().getNBins(); bin++) {
                
                H1F h = h2.sliceX(bin);
                h.setTitleX("distance (mm)");
                h.setTitleY("count");

                double mean = h.getMean();
                double width = h.getRMS();
                GraphErrors func = null;
                Pair<CrystalBall, GraphErrors> fitResult = fit_with_crystal_ball(h, "time2distance, wire " + (i+1) + "/576, slice " + bin + " in analyse_wire_histograms()...");
                if (fitResult != null) {
                    CrystalBall cb = fitResult.getFirst();
                    // System.out.println("* Fit result : time2distance layer + " + i + ",  bin " + bin);
                    // cb.print();
                    func = fitResult.getSecond();
                    mean = cb.getMu();
                    width = cb.getSigma();
                }

                // store results
                //if (width < 0.5){
                if (true){
                    gr.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, width);
                    gr_bis.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, 0);
                }
                h.setTitle(String.format("mean : %f , width %f", mean, width));

                // output
                save_histo1D(h, func, String.format("%s/time2distance-wire-L%dW%d-slice-%d.pdf", slicesDir, identifier.layer, identifier.component, bin));
            } // end loop over bin

            // Data point to be fit
            ArrayList<Double> xList = new ArrayList<>();
            ArrayList<Double> yList = new ArrayList<>();
            // xList.add(0.0);
            // yList.add(0.0);
            for (int pt = 0; pt < gr_bis.getDataSize(0)-2; pt++) { // exclude the 2 end points
                double x = gr_bis.getDataX(pt);
                double y = gr_bis.getDataY(pt);
                if (x < 230) { // filtering
                    xList.add(x);
                    yList.add(y);
                }
            }
            double[] xData = new double[xList.size()];
            double[] yData = new double[xList.size()];

            for (int bin = 0; bin < xData.length; bin++) {
                xData[bin] = xList.get(bin);
                yData[bin] = yList.get(bin);
            }

            // Fit
            GraphErrors func = null;
            double[] initialValues = results.time2distance.get(i); // use the actual t2d values to initialise the fit

            try {
                double[] fittedParameters = fit_time2distance(xData, yData, initialValues.clone());
                int Npts = 1000;
                double tmin = 0;
                double tmax = 250;
                func = new GraphErrors();
                func.setLineColor(2);
                func.setLineThickness(1);
                func.setMarkerSize(0);
                for (int bin = 0; bin < Npts; bin++) {
                    double time = tmin + (tmax-tmin)*bin/(Npts-1);
                    func.addPoint(time, eval_t2d(time, fittedParameters), 0, 0);
                }
                t2d_params_list.add(fittedParameters);
            } catch (Exception e) {
                System.out.println("Fail to fit time2disatnce: wire L" + identifier.layer + "W" + identifier.component);
                // System.out.println("Error: " + e.getMessage());
                // e.printStackTrace();
                t2d_params_list.add(initialValues.clone());
            }
            
            int canvas_num = i / 12; // between 0 and 47
            int canvas_frame = i % 12; // between 0 and 11
            canvas2.get(canvas_num).cd(canvas_frame);
            canvas2.get(canvas_num).getPad(canvas_frame).setPalette("kBird");
            canvas2.get(canvas_num).getPad(canvas_frame).getAxisFrame().setDrawAxisZ(false);
            h2_initial.setTitle("time2distance");
            canvas2.get(canvas_num).draw(h2_initial);
            if (func != null) canvas2.get(canvas_num).draw(func, "same L");
            canvas2.get(canvas_num).draw(gr_bis, "same P");
            canvas2.get(canvas_num).getPad(canvas_frame).getAxisY().setRange(0, 3);
            canvas2.get(canvas_num).getPad(canvas_frame).getAxisX().setRange(0, 250);
            canvas2.get(canvas_num).getPad(canvas_frame).getAxisFrame().setDrawAxisZ(false); // again

        }
        // save histograms
        for (int i = 0; i < 48; i++) {
            canvas2.get(i).save(wireOutDir + String.format("/time2distance_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
            System.out.println(String.format("time2distance_iter_%d_w%d-%d.pdf", niter, 12*i, 12*(i+1)-1));
        }
        results.time2distance = t2d_params_list;

    }


    static void fill_histos(Histos histos, DataEvent event, AlertTrackSelector analyser, ALERTEngine alertEngine, AHDCEngine ahdcEngine, boolean flag_do_fit) {
        //if (analyser.hasElasticElectron(event)) {
        if (analyser.hasGoodTrack(event)) {
            // // Retrieve the kinematics of the electron
            // ParticleRow electron = analyser.getElectron();

            // Run engine to fill the residual_LR
            ahdcEngine.processDataEventUser(event);
            if (flag_do_fit) {
                alertEngine.processDataEventUser(event);
            } else {
                //alertEngine.processDataEventProjOnly(event, AHDCdet);
                
            }
            

            // Data analysis
            DataBank trackBank = event.getBank("AHDC::kftrack");
            DataBank hitBank = event.getBank("AHDC::hits");
            DataBank recBank = event.getBank("REC::Particle");
            //recBank.show();

            // Look at residuals in elastic tracks
             ArrayList<int[]> trackRows = analyser.getAhdcKFTrackRows();
            
            for (int i = 0; i < trackRows.size(); i++) {
                int[] rows = trackRows.get(i);
                int track_row = rows[1];
                int trackid = trackBank.getInt("trackid", track_row);
                
                double track_vz = trackBank.getFloat("z" , track_row); // mm
                double px = trackBank.getFloat("px", track_row); // MeV
                double py = trackBank.getFloat("py", track_row); // MeV
                double pz = trackBank.getFloat("pz", track_row); // MeV
                double p = Math.sqrt(Math.pow(px, 2)+Math.pow(py, 2)+Math.pow(pz, 2));
                double track_theta_deg = Math.toDegrees(Math.acos(pz/p));
                
                // Layer or wire histograms
                for (int j = 0; j < hitBank.rows(); j++) {
                    if (trackid == hitBank.getInt("trackid", j)) {
                        int layer = 10*hitBank.getByte("superlayer", j) + hitBank.getByte("layer", j);
                        double residual = hitBank.getDouble("residual", j); // mm
                        double residual_LR = hitBank.getDouble("residual_LR", j); // mm, residual_LR stocked here for dev purpose
                        int component = hitBank.getInt("wire", j);
                        AhdcWireId wireId = new AhdcWireId(1, layer, component);
                        double time = hitBank.getDouble("time", j);
                        double doca = hitBank.getDouble("doca", j);
                        double distance = doca - residual;

                        double vz_hit = track_vz + AhdcWireId.layer2Radius(layer)/Math.tan(Math.toRadians(track_theta_deg));

                        // Fill histograms per layers
                        histos.h1_residual_LR_per_layers.get(AhdcWireId.layer2number(layer)).fill(residual_LR);
                        histos.h1_residual_per_layers.get(AhdcWireId.layer2number(layer)).fill(residual);
                        histos.h2_corr_vz_residual_LR_per_layers.get(AhdcWireId.layer2number(layer)).fill(vz_hit, residual_LR);
                        histos.h2_time2distance_per_layers.get(AhdcWireId.layer2number(layer)).fill(time, distance);
                            // integrated
                        histos.h1_residual_LR_per_layers.get(0).fill(residual_LR);
                        histos.h1_residual_per_layers.get(0).fill(residual);
                        histos.h2_corr_vz_residual_LR_per_layers.get(0).fill(vz_hit, residual_LR);
                        histos.h2_time2distance_per_layers.get(0).fill(time, distance);

                        // Fill histos per wires
                        histos.h1_residual_LR_per_wires.get(wireId.num).fill(residual_LR);
                        histos.h2_corr_vz_residual_LR_per_wires.get(wireId.num).fill(vz_hit, residual_LR);
                        histos.h2_time2distance_per_wires.get(wireId.num).fill(time, distance);

                    } // end select hits for this trackid
                }  // end loop over hits

                // Integrated distributions
                histos.h1_track_theta.fill(track_theta_deg);

                // link with the electron
                int electron_row = rows[0];
                    
                    // delta phi
                double phi_track_rad = Math.atan2(trackBank.getFloat("py", track_row), trackBank.getFloat("px", track_row));
                double phi_electron_rad = Math.atan2(recBank.getFloat("py", electron_row), recBank.getFloat("px", electron_row));
                double delta_phi_rad = diffAngles(phi_electron_rad, phi_track_rad);
                delta_phi_rad = diffAngles(delta_phi_rad, Math.PI); // remove 180°
                histos.h1_track_delta_phi.fill(Math.toDegrees(delta_phi_rad));
                    // delta vz
                double delta_vz = recBank.getFloat("vz", electron_row) - 0.1*track_vz; // cm
                histos.h1_delta_vz.fill(delta_vz);

            } // end loop over good track rows

        } // end event has good track
    }

    static ArrayList<String> verify_files(ArrayList<String> inFiles) {
        ArrayList<String> goodfiles = new ArrayList<>();
        for (String file : inFiles) {
            if (Files.exists(Path.of(file))) {
                goodfiles.add(file);
            }
        }
        return goodfiles;
    }

    static void check_output_dir(String outDir) {
        if (!Files.exists(Path.of(outDir))) {
            try {
                Files.createDirectories(Path.of(outDir));
                System.out.println("Output directory created and set to : " + outDir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     *  Here all layers are rotated. Equivalent to {@link #computeNewLayerAngles(double[], double[], int)} for nlayers = 8
     */
    static void computeNewLayerAngles(double[] angles, double[] residuals) {
        for (int i = 0; i < angles.length; i++) {
            double alphaRad = residuals[i]/AhdcWireId.layerNum2Radius(i+1);
            angles[i] = angles[i] - 0.5*Math.toDegrees(alphaRad);
        }
    }

    /**
     * Rotate the layers with the worst shifts
     * @param angles angles to be updated 
     * @param residuals residuals to estimate the new correction angles
     * @param nlayers number of layers to be rotated.
     */
    static void computeNewLayerAngles(double[] angles, double[] residuals, int nlayers) {
        double[] abs_residuals = new double[residuals.length];
        for (int i = 0; i < residuals.length; i ++) {
            abs_residuals[i] = Math.abs(residuals[i]);
        }
        /// --- sort values according to the residuals
        int[] sortedIndices = sortIndices(abs_residuals, false);
        /// --- only look at the nlayers worst residuals
        for (int loop = 0; loop < sortedIndices.length && loop < nlayers; loop++) {
            int i = sortedIndices[loop];
            double alphaRad = residuals[i]/AhdcWireId.layerNum2Radius(i+1);
            angles[i] = angles[i] - 0.5*Math.toDegrees(alphaRad);
        }
    }

    /**
     * As for {@link #computeNewLayerAngles(double[], double[], int) but for wires}
     * @param angles angles to be updated 
     * @param residuals residuals to estimate the new correction angles
     * @param nwires number of wires to be rotated.
     */
    static void computeNewWireAngles(double[] angles, double[] residuals, int nwires) {
        double[] abs_residuals = new double[residuals.length];
        for (int i = 0; i < residuals.length; i ++) {
            abs_residuals[i] = Math.abs(residuals[i]);
        }
        /// --- sort values according to the residuals
        int[] sortedIndices = sortIndices(abs_residuals, false);
        /// --- only look at the nlayers worst residuals
        for (int loop = 0; loop < sortedIndices.length && loop < nwires; loop++) {
            int i = sortedIndices[loop];
            AhdcWireId identifier = new AhdcWireId(i);
            double radius = AhdcWireId.layer2Radius(identifier.layer);
            double alphaRad = residuals[i]/radius;
            angles[i] = angles[i] - 0.5*Math.toDegrees(alphaRad);
        }
    }

    static void printRotationAngles(ResultsOverIterations results, boolean flag_wire) {
        if (!flag_wire) {
            double[] start = results.layer_angles_start;
            double[] end = results.layer_angles_end;
            System.out.println("----------------------------------------------");
            System.out.println("Correction angles to be applied");
                System.out.printf("   %5s   %6s    %6s\n", "layer", "start", "end");
            for (int i = 0; i < start.length; i++) {
                System.out.printf("   %5d   %6.4f    %6.4f\n", AhdcWireId.number2layer(i+1), start[i], end[i]);
            }
            System.out.println("----------------------------------------------");
        } else {
            System.out.println("to be done later...");
        }
    }

    /**
     * New method to update rotation angles
     * @param results
     */
    static void computeNewWireAngles(ResultsOverIterations results) {
        double[] angles = results.wire_angles;
        double[] residuals = results.wire_residuals;
        for (int i = 0; i < angles.length; i++) {
            AhdcWireId wireId = new AhdcWireId(i);
            double alphaRad = residuals[i]/AhdcWireId.layer2Radius(wireId.layer);
            angles[i] = angles[i] - 0.5*Math.toDegrees(alphaRad);
        }
    }

    /**
     * 
     * @param vector data to be sorted
     * @return the array of id in the right order
     */
    public static int[] sortIndices(double[] vector, boolean rising) {
        class Data implements Comparable<Data> {
            private int id;
            private double value;
            public Data(int _id, double _value) { id = _id; value = _value;}
            public int compareTo(Data data) {
                return Double.compare(this.value, data.value);
            }
            public int getInt() { return id;}
            //public double getValue() {return value;}
        }

        ArrayList<Data> list = new ArrayList<>();
        for (int i = 0; i < vector.length; i++) {
            list.add(new Data(i, vector[i]));
        }
        if (rising) { // rising
            list.sort(null);
        } else { // inverse order
            list.sort(Collections.reverseOrder());
        }

        int[] indices = new int[vector.length];
        for (int i = 0; i < list.size(); i++) {
            indices[i] = list.get(i).getInt();
        }

        return indices;
    }

  /**
   * Generate the AHDC geometry with specific correction angles
   * @param _wire_angles_start rotation to be applied to the start of the AHDC wires
   * @param _wire_angles_end rotation to be applied to the end of the AHDC wires
   * @return AHDC geometry
   */
    // static AlertDCDetector generateAhdcGeometry(double[] _wire_angles_start, double[] _wire_angles_end) {
    //     AlertDCFactory factory = new AlertDCFactory();
    //     factory.setWireCorrectionAngles(_wire_angles_start, _wire_angles_end);
    //     return factory.createDetectorCLAS(new DatabaseConstantProvider());
    // }

    /**
     * Return the one by one sum of the component sof a and b
     * @param a array
     * @param b array
     * @return new array
     */
    static public double[] sum_array(double[] a, double[] b) {
        if (a.length != b.length) {
            return null;
        } else {
            double[] res = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                res[i] = a[i] + b[i];
            }
            return res;
        }
    }

    /**
     * Assign the value a to all components of the array v
     * @param v array
     * @param a double
     */
    static public void assign_value(double[] v, double a) {
        if (v == null) return ;
        for (int i = 0; i < v.length; i++) {
            v[i] = a;
        }
    }

    /**
     * Inverse operation of {@link #doLayerRotations(AlertDCDetector, double[])}
     * @param AHDCdet AHDC detector to be rotated
     * @param layer_angles rotation angles for the layers
     */
    static void undoLayerRotations(AlertDCDetector AHDCdet, double[] layer_angles) {
        for (int num = 0; num < 8; num++) {
            // rotate AHDCdet
            int layer = AhdcWireId.number2layer(num+1);
            int sl = layer / 10;
            int l = layer % 10;
            for (int i = 0; i < AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getNumComponents(); i++) {
                AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getComponent(i+1); // numbering starts at 1
                wire.rotateZ(Math.toRadians(-layer_angles[num]));
            }
        }
    }

    /**
     * Rotate the layers of the AHDC
     * @param AHDCdet AHDC detector to be rotated
     * @param layer_angles rotation angles for the layers
     */
    static void doLayerRotations(AlertDCDetector AHDCdet, double[] layer_angles) {
        System.out.println("\033[1;32m > Rotate AHDC detector \033[0m");
        for (int num = 0; num < 8; num++) {
            // rotate AHDCdet
            int layer = AhdcWireId.number2layer(num+1);
            int sl = layer / 10;
            int l = layer % 10;
            for (int i = 0; i < AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getNumComponents(); i++) {
                AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getComponent(i+1); // numbering starts at 1
                wire.rotateZ(Math.toRadians(layer_angles[num]));
            }
            //System.out.println("   layer " + layer + " --> will be rotated by " + layer_angles[num] + " deg");
            System.out.printf("   layer " + layer + " --> will be rotated by %.2f deg \n", layer_angles[num]);
        }
    }

    /**
     * Layer alignment analysis
     * @param inFiles HIPO files
     * @param outDir output directory
     * @param flag_do_fit does KF perform a fit ?
     * @throws IOException 
     * @throws InterruptedException 
     */
    static void layer_alignment(ArrayList<String> inFiles, String outDir, boolean flag_do_fit, Supplier<AlertTrackSelector> trackSelectorFactory, HistoAnalyser... analysers) throws IOException, InterruptedException {

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        /// --- Define initial geometry according the values in ResultsOverIterations
        //AlertDCDetector AHDCdet = generateAhdcGeometry(layerAngles2WireAngles(results.layer_angles_start), layerAngles2WireAngles(results.layer_angles_end));

        check_output_dir(outDir + "/initialConfig");
        results.screenshot(outDir + "/initialConfig");

        /// --- Global observables
        GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
        g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
        g0.setTitleX("iterations");
        g0.setTitleY("cost");

        ArrayList<GraphErrors> gr_slopes = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            GraphErrors gr = new GraphErrors();
            //gr.setMarkerSize(0);
            gr.setMarkerColor(i+1);
            gr.setLineColor(i+1);
            gr.setTitleX("iterations");
            gr.setTitleY("layer slope");
            gr_slopes.add(gr);
        }

        ArrayList<GraphErrors> gr_constants = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            GraphErrors gr = new GraphErrors();
            //gr.setMarkerSize(0);
            gr.setMarkerColor(i+1);
            gr.setLineColor(i+1);
            gr.setTitleX("iterations");
            gr.setTitleY("layer constant");
            gr_constants.add(gr);
        }

        /// --- Loop over criteria
        int niter = 0;
        double value = 1e10;
        while (niter < 25) {
            niter++;
            // run iteration
            System.out.println("\033[1;32m ################################ \033[0m");
            System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
            System.out.println("\033[1;32m ################################ \033[0m");

            run(niter, results, inFiles, outDir, +75, flag_do_fit, trackSelectorFactory, analysers);

            // running criteria
            value = 0;
            int N = 0;
            for (int i = 0; i < results.layer_residuals.length; i++) {
                value += Math.pow(results.layer_residuals[i], 2);
                N++;
            }
            value = Math.sqrt(value) / N;
            g0.addPoint(niter, value, 0, 0);
            System.out.println("\033[1;31m =======> convergence criteria : " + value + "\033[0m");

            /// --- Update rotation angles
            computeNewLayerAngles(results.layer_angles_start, results.layer_residuals_start, 8); // wire start 
            computeNewLayerAngles(results.layer_angles_end, results.layer_residuals_end, 8); // wire end
            printRotationAngles(results, false);

            ///--- output
            for (int i = 0; i < 8; i++) {
                gr_slopes.get(i).addPoint(niter, results.layer_residuals_slope[i], 0, 0);
                gr_constants.get(i).addPoint(niter, results.layer_residuals_constant[i], 0, 0);
            }

            /// --- save ResultOverIterations and update ccdb table for the new geometry
            results.update_ccdb(outDir + "/layers/iter-" + niter, CCDB_TYPE.LAYER);

        } // end loop over criteria / nb iterations
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
        System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

        EmbeddedCanvas c1 = new EmbeddedCanvas(1200, 900);
        EmbeddedCanvas c2 = new EmbeddedCanvas(1200, 900);
        for (int i = 0; i < 8; i++) {
            if (i == 0) {
                c1.draw(gr_slopes.get(i), "P");
                c2.draw(gr_constants.get(i), "P");
            } else {
                c1.draw(gr_slopes.get(i), "same P");
                c2.draw(gr_constants.get(i), "same P");
            }
        }

        c1.save(outDir + "/summary-layer-slopes-over-iterations.pdf");
        c2.save(outDir + "/summary-layer-constants-over-iterations.pdf");

    }

    /**
     * Time2distance calibration
     * @param inFiles HIPO files
     * @param outDir output directory
     * @param flag_do_fit does KF perform a fit ?
     * @throws IOException 
     * @throws InterruptedException 
     */
    static void time2distance_calibration(ArrayList<String> inFiles, String outDir, boolean flag_do_fit, Supplier<AlertTrackSelector> trackSelectorFactory, HistoAnalyser... analysers) throws IOException, InterruptedException {

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        check_output_dir(outDir + "/initialConfig");
        results.screenshot(outDir + "/initialConfig");

        /// --- Global observables
        GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
        g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
        g0.setTitleX("iterations");
        g0.setTitleY("cost");

        ArrayList<GraphErrors> gr_slopes = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            GraphErrors gr = new GraphErrors();
            //gr.setMarkerSize(0);
            gr.setMarkerColor(i+1);
            gr.setLineColor(i+1);
            gr.setTitleX("iterations");
            gr.setTitleY("layer slope");
            gr_slopes.add(gr);
        }

        ArrayList<GraphErrors> gr_constants = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            GraphErrors gr = new GraphErrors();
            //gr.setMarkerSize(0);
            gr.setMarkerColor(i+1);
            gr.setLineColor(i+1);
            gr.setTitleX("iterations");
            gr.setTitleY("layer constant");
            gr_constants.add(gr);
        }

        /// --- Loop over criteria
        int niter = 0;
        double value = 1e10;
        while (niter < 25) {
            niter++;
            // run iteration
            System.out.println("\033[1;32m ################################ \033[0m");
            System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
            System.out.println("\033[1;32m ################################ \033[0m");

            run(niter, results, inFiles, outDir, +75, flag_do_fit, trackSelectorFactory, analysers);

            // running criteria
            value = 0;
            int N = 0;
            for (int i = 0; i < results.layer_residuals.length; i++) {
                value += Math.pow(results.layer_residuals[i], 2);
                N++;
            }
            value = Math.sqrt(value) / N;
            g0.addPoint(niter, value, 0, 0);
            System.out.println("\033[1;31m =======> convergence criteria : " + value + "\033[0m");

            ///--- output
            for (int i = 0; i < 8; i++) {
                gr_slopes.get(i).addPoint(niter, results.layer_residuals_slope[i], 0, 0);
                gr_constants.get(i).addPoint(niter, results.layer_residuals_constant[i], 0, 0);
            }

            /// --- save ResultOverIterations and update ccdb table for the new geometry
            results.update_ccdb(outDir + "/layers/iter-" + niter, CCDB_TYPE.T2D);

        } // end loop over criteria / nb iterations
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
        System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

        EmbeddedCanvas c1 = new EmbeddedCanvas(1200, 900);
        EmbeddedCanvas c2 = new EmbeddedCanvas(1200, 900);
        for (int i = 0; i < 8; i++) {
            if (i == 0) {
                c1.draw(gr_slopes.get(i), "P");
                c2.draw(gr_constants.get(i), "P");
            } else {
                c1.draw(gr_slopes.get(i), "same P");
                c2.draw(gr_constants.get(i), "same P");
            }
        }

        c1.save(outDir + "/summary-layer-slopes-over-iterations.pdf");
        c2.save(outDir + "/summary-layer-constants-over-iterations.pdf");

    }

    /**
     * Code obselète: à modifier pour épouser la nouvelle philosophie.
     *                  --> modifier la ccdb et créer la géométrie à partir de la ccdb
     *                  --> cf. layer_alignment() for example
     * 
     * To be done later.
     * 
     * Routine to scan ahdc position with respect to CLAS
     * @param args
     * @param flag_do_fit if true, a KF fit is performed, if false, a simple propagation is performed
     */
    static void scan_ahdc_position(String[] args, boolean flag_do_fit) {

        /// --- Load inputs from options
        fOptions options = new fOptions("-i", "-o");
        options.LoadOptions(args);
        options.Show();

        /// --- Verify input files are not empty
        ArrayList<String> inFiles = verify_files(options.GetValues("-i"));
        if (inFiles.size() == 0) {
            System.out.println("Please provide inputs files using the option: -i");
            return;
        }

        /// --- Check that the output dir exists or create a new one
        String outDir0 = options.GetValue("-o");
        check_output_dir(outDir0);

        /// --- Obsevables
        GraphErrors g_mean = new GraphErrors("mean-angle-over-position");
        g_mean.setTitle("Mean angle over clas alignment");
        g_mean.setTitleX("clas alignment");
        g_mean.setTitleY("mean angle");
        GraphErrors g_dev = new GraphErrors("angle-deviation-over-position");
        g_dev.setTitle("Angle deviation over clas alignment");
        g_dev.setTitleX("clas alignment");
        g_dev.setTitleY("angle deviation");
        GraphErrors g_mean_vs_dev = new GraphErrors("mean-angle-vs-deviation");
        g_mean_vs_dev.setTitle("Mean vs deviation");
        g_mean_vs_dev.setTitleX("mean angle");
        g_mean_vs_dev.setTitleY("angle deviation");

        /// --- Loop over CLAS position with respect to CLAS
        /// Scan from 30 to 100 mm
        for (int step = 52; step < 100; step++) {

            String outDir = outDir0 + "/" + step;
            check_output_dir(outDir);

            double clas_alignment = 1.0*step; // mm

            System.out.println("\033[1;33m >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \033[0m");
            System.out.println("\033[1;33m # CLAS alignment : " + clas_alignment + "\033[0m");
            System.out.println("\033[1;33m <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\033[0m");

            // --- Results over iterations
            ResultsOverIterations results = new ResultsOverIterations();

            /// --- Define initial geometry parameters
            AlertDCDetector AHDCdet = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

            /// --- rotate AHDC detector
            System.out.println(" Initial rotation angles : AHDC detector");
            doLayerRotations(AHDCdet, results.layer_angles);

            /// --- Global observables
            GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
            g0.setTitle("Convergence, clas_alignment : " + clas_alignment);
            g0.setTitleX("iterations");
            g0.setTitleY("convergence criteria");
                
            /// --- Loop over criteria
            int niter = 0;
            double value = 1e10;
            //while (niter < 12) {
            while (niter < 25) {
                niter++;
                // run iteration
                System.out.println("\033[1;32m ########### Start iteration : " + niter + "\033[0m clas_alignment : " + clas_alignment);

                //run(niter, results, inFiles, outDir, AHDCdet, clas_alignment, false, flag_do_fit);
                // we will not be able to run the code in this condition
                // need to adapt the doLayerRotations to rely on a sqlite file
                // let fix layer and wire alignmpent first

                // running criteria
                value = 0;
                int N = 0;
                for (int i = 0; i < results.layer_residuals.length; i++) {
                    value += Math.pow(results.layer_residuals[i], 2);
                    N++;
                }
                value = Math.sqrt(value) / N;
                g0.addPoint(niter, value, 0, 0);
                System.out.println("\033[1;31m =======> convergence criteria : " + value + "\033[0m");

                // look at the distribution of angles
                if (value < 2*1e-3 || niter >= 25) {
                    double sum = 0;
                    double sum2 = 0;
                    int npts = 0;
                    for (int i = 0; i < results.layer_angles.length; i++) {
                        npts++;
                        sum += results.layer_angles[i];
                        sum2 += Math.pow(results.layer_angles[i], 2);
                    }
                    double mean = sum/npts;
                    double var = (sum2/npts) - Math.pow(mean, 2);
                    double deviation = Math.sqrt(var);
                    g_mean.addPoint(clas_alignment, mean, 0, 0);
                    g_dev.addPoint(clas_alignment, deviation, 0, 0);
                    g_mean_vs_dev.addPoint(mean, deviation, 0, 0);

                    System.out.println("\033[1;33m ************ angle mean : " + mean + " , deviation : " + deviation + "\033[0m , clas_alignment : " + clas_alignment);
                    break;
                }
                
                /// --- Undo AHDC rotation before applying new rotation angles to prevent accumulation
                undoLayerRotations(AHDCdet, results.layer_angles);

                /// --- Update rotation angles
                computeNewLayerAngles(results.layer_angles, results.layer_residuals);

                /// --- Rotate AHDC detector
                doLayerRotations(AHDCdet, results.layer_angles);
                
            } // end loop over criteria / nb iterations
            EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
            c0.draw(g0);
            c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
            System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

        } // end loop over steps      
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g_mean);
        c0.save(outDir0 + "/summary-mean-angle.pdf");
        System.out.println(outDir0 + "/summary-mean-angle.pdf");
        EmbeddedCanvas c1 = new EmbeddedCanvas(1200, 900);
        c1.draw(g_dev);
        c0.save(outDir0 + "/summary-angle-deviation.pdf");
        System.out.println(outDir0 + "/summary-angle-deviation.pdf");
        EmbeddedCanvas c2 = new EmbeddedCanvas(1200, 900);
        c2.draw(g_mean_vs_dev);
        c2.save(outDir0 + "/summary-mean-vs-dev.pdf");
        System.out.println(outDir0 + "/summary-mean-vs-dev.pdf");

    } // end scan ahdc position

    /**
     * New AHDC position scan: the oversable is delta vz versus clas alignment
     * @param inFiles HIPO files
     * @param outDir output directory
     * @param flag_do_fit does KF perform a fit ?
     * @throws IOException 
     * @throws InterruptedException 
     */
    static void new_ahdc_position_scan(ArrayList<String> inFiles, String outDir, boolean flag_do_fit, Supplier<AlertTrackSelector> trackSelectorFactory, HistoAnalyser... analysers) throws IOException, InterruptedException {

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        check_output_dir(outDir + "/initialConfig");
        results.screenshot(outDir + "/initialConfig");

        /// --- Obsevables
        GraphErrors g_mean = new GraphErrors("mean-angle-over-position");
        g_mean.setTitle("Mean delta vz versus clas alignment");
        g_mean.setTitleX("clas alignment");
        g_mean.setTitleY("mean delta vz");
        GraphErrors g_width = new GraphErrors("angle-deviation-over-position");
        g_width.setTitle("Width delta vz versus clas alignment");
        g_width.setTitleX("clas alignment");
        g_width.setTitleY("width delta vz");

        BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "/results.txt"));
        writer.write("# clas alignment, mean delta vz, width delta vz");
        writer.newLine();
        writer.newLine();

        /// --- Loop over criteria
        for (int step = 30; step < 80; step++) {

            String stepDir = outDir + "/" + step;
            check_output_dir(stepDir);

            double clas_alignment = 1.0*step; // mm
            results.clas_alignment = clas_alignment;

            // run iteration
            System.out.println("\033[1;32m ########### Start iteration : \033[0m clas_alignment : " + clas_alignment);

            run(0, results, inFiles, stepDir, clas_alignment, flag_do_fit, trackSelectorFactory, analysers);

            double mean = results.mean_delta_vz;
            double width = results.width_delta_vz;
            
            g_mean.addPoint(clas_alignment, mean, 0,0);
            g_width.addPoint(clas_alignment, width, 0,0);

            String line = String.format("%f   %f   %f", clas_alignment, mean, width);
            writer.write(line);
            writer.newLine();

        } // end loop over criteria / nb iterations
        writer.close();

        { // Mean
            EmbeddedCanvas c = new EmbeddedCanvas(800, 600);
            c.draw(g_mean);
            c.save(outDir + "/delta_vz_mean_versus_clas_alignment.pdf");
        }
        { // Width
            EmbeddedCanvas c = new EmbeddedCanvas(800, 600);
            c.draw(g_width);
            c.save(outDir + "/delta_vz_width_versus_clas_alignment.pdf");
        }
    }

    /**
     * Wire alignment analysis
     * @param inFiles HIPO files
     * @param outDir output directory
     * @param flag_do_fit does KF perform a fit ?
     * @throws IOException 
     * @throws InterruptedException 
     */
    static void wire_alignment(ArrayList<String> inFiles, String outDir, boolean flag_do_fit, Supplier<AlertTrackSelector> trackSelectorFactory, HistoAnalyser... analysers) throws IOException, InterruptedException {

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        /// --- Define initial geometry according the values in ResultsOverIterations (start from layer alignment)
        // results.wire_angles_start = layerAngles2WireAngles(results.layer_angles_start); // not to be changed, based on the layer alignment
        // results.wire_angles_end   = layerAngles2WireAngles(results.layer_angles_end); // not to be changed, based on the layer alignment
        // //AlertDCDetector AHDCdet = generateAhdcGeometry(results.wire_angles_start, results.wire_angles_end);
        // //assign_value(results.wire_angles, 0.0); // is defined with respect to the layer alignment
        // results.wire_angles = init_wire_angles(); // read last alignment from values save in wire_angles.txt
        // assign_value(results.wire_residuals, 0.0);

        check_output_dir(outDir + "/initialConfig");
        results.screenshot(outDir + "/initialConfig");

        /// --- Global observables
        GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
        g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
        g0.setTitleX("iterations");
        g0.setTitleY("cost");

        /// --- Loop over criteria
        int niter = 0;
        double value = 1e10;
        //while (niter < 12) {
        //while (value > 1*1e-3 && niter < 25) {
        while (niter < 25) {
            niter++;
            // run iteration
            System.out.println("\033[1;32m ################################ \033[0m");
            System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
            System.out.println("\033[1;32m ################################ \033[0m");

            run(niter, results, inFiles, outDir, +75, flag_do_fit, trackSelectorFactory, analysers);

            // running criteria
            value = 0;
            int N = 0;
            for (int i = 0; i < results.wire_residuals.length; i++) {
                value += Math.pow(results.wire_residuals[i], 2);
                N++;
            }
            value = Math.sqrt(value) / N;
            g0.addPoint(niter, value, 0, 0);
            System.out.println("\033[1;31m =======> convergence criteria : " + value + "\033[0m");

            /// --- Update rotation angles
            // computeNewWireAngles(results.wire_angles_start, results.wire_residuals_start, 576); // wire start 
            // computeNewWireAngles(results.wire_angles_end, results.wire_residuals_end, 576); // wire end
            computeNewWireAngles(results.wire_angles, results.wire_residuals, 576);
            //printRotationAngles(results, false);

            /// --- Update AHDC geometry
            //AHDCdet = generateAhdcGeometry(results.wire_angles_start, results.wire_angles_end);
            //AHDCdet = generateAhdcGeometry(sum_array(results.wire_angles_start, results.wire_angles), sum_array(results.wire_angles_end, results.wire_angles));

            /// --- save ResultOverIterations
            results.update_ccdb(outDir + "/wires/iter-" + niter, CCDB_TYPE.WIRE);

        } // end loop over criteria / nb iterations
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
        System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

    }

    /**
     * Save histogram as pdf
     * @param h histogram
     * @param gr fit function
     * @param filename output name (full path)
     */
    public static void save_histo1D(H1F h, GraphErrors gr, String filename) {
        EmbeddedCanvas c = new EmbeddedCanvas(800, 600);
        c.draw(h);
        if (gr != null)
            c.draw(gr, "same L");
        c.save(filename);
    }

    /**
     * Fit the hsitogram with a cristal ball distribution. Return null if the fit failed
     * @param h histogram to be fitted
     * @return
     */
    public static Pair<CrystalBall, GraphErrors> fit_with_crystal_ball(H1F h, String debug_name) {
        // initialize fit parameters
        float[] data = h.getData();
        double amp = h.getMax();
        int binMax = h.getMaximumBin();
        int binHalfLeft = -1;
        int binHalfRight = -1;
        int binTenthLeft = -1;
        int binTenthRight = -1;
        for (int bin = 0; bin < data.length-1; bin++) {
            if (bin < binMax) {
                if (data[bin] < 0.5*amp && data[bin+1] > 0.5*amp) {
                    binHalfLeft = bin;
                }
                if (data[bin] < 0.1*amp && data[bin+1] > 0.1*amp) {
                    binTenthLeft = bin;
                }
            } else {
                if (data[bin] > 0.5*amp && data[bin+1] < 0.5*amp) {
                    binHalfRight = bin;
                }
                if (data[bin] > 0.1*amp && data[bin+1] < 0.1*amp) {
                    binTenthRight = bin;
                }
            }
        }
        double x1 = h.getxAxis().getBinCenter(binHalfLeft); // value at half amplitude
        double x2 = h.getxAxis().getBinCenter(binHalfRight); // value at half amplitude
        double xpeak = h.getxAxis().getBinCenter(h.getMaximumBin());
        double x01 = h.getxAxis().getBinCenter(binTenthLeft); // value at 0.1* amplitude
        double x02 = h.getxAxis().getBinCenter(binTenthRight); // value at 0.1* amplitude

        double sigma = (x2-x1)/2;
        double xmin = xpeak - 1.4*sigma;
        double xmax = xpeak + 1.4*sigma;
        
        double alpha0 = 0;
        double leftWidth  = xpeak - x01;
        double rightWidth = x02 - xpeak;
        CrystalBall.QueueSide cbSide = CrystalBall.QueueSide.LEFT; // per defaulf, queue à gauche
        if (leftWidth > rightWidth) { // right asymmetry
            alpha0 = leftWidth/sigma;
        } else {
            alpha0 = rightWidth/sigma;
            cbSide = CrystalBall.QueueSide.RIGHT;
        }

        // perform fit
        CrystallBallFitter cbFitter = new CrystallBallFitter();
        //cbFitter.setAlphaParameter(Math.abs(alpha0), Math.max(0.7*Math.abs(alpha0), 1e-3), + 1.3*Math.abs(alpha0));
        cbFitter.setAlphaParameter(alpha0, 0.3, 5);
        cbFitter.setNpowerParameter(5, 3, 100); // 3 is a standard value
        cbFitter.setMuParameter(xpeak, xmin, xmax);
        cbFitter.setSigmaParameter(sigma, Math.max(0.6*sigma, 1e-6), 1.1*sigma);
        cbFitter.setAmplitudeParameter(amp, 0.8*amp, 1.2*amp);
        cbFitter.setQueueSide(cbSide);
        //cbFitter.print();

        try {
            Pair<CrystalBall, GraphErrors> fitResult = cbFitter.fit(h, xmin, xmax);
            return fitResult;
        } catch (Exception e) {
            System.out.println("Crystal Ball fit failed for histogram: "+ debug_name);
            // cbFitter.print();
            // System.out.println("Error: " + e.getMessage());
            // e.printStackTrace();
            return null;
        }
    }

    /**
     * Current time2distance parametrization of Michael. Ref. coatjava
     * @param time variable
     * @param params parameters
     * @return the distance with respect to the AHDC wire
     */
    static double eval_t2d(double time, double[] params) {
        double p1_int = params[0];
        double p1_slope = params[1];
        double p2_int = params[2];
        double p2_slope = params[3];
        double p3_int = params[4];
        double p3_slope = params[5];
        double t1_x0 = params[6];
        double t1_width = params[7];
        double t2_x0 = params[8];
        double t2_width = params[9];

        double p1 = p1_int + p1_slope*time;
        double p2 = p2_int + p2_slope*time;
        double p3 = p3_int + p3_slope*time;

        double t1 = 1.0/(1.0 + Math.exp(-(time - t1_x0)/t1_width));
        double t2 = 1.0/(1.0 + Math.exp(-(time - t2_x0)/t2_width));

        //return p1*(1.0 - t1) + t1*(1.0 - t2)*p2 + t1*t2*p3;
        return p1*(1.0 - t1) + t1*(1.0 - t2)*p2 + t2*p3;
    }

    /**
     * Return the parameters of the time2distance used by by {@link #eval_t2d(double, double[])}
     * @param xData
     * @param yData
     * @param initialValues
     * @return
     */
    public static double[] fit_time2distance(double[] xData, double[] yData, double[] initialValues) {

        // --- Compute a new vector yDataModelled acoording to the model
        MultivariateJacobianFunction model = new MultivariateJacobianFunction() {
            
            public Pair<RealVector, RealMatrix> value(final RealVector parameters) {

                double[] params = parameters.toArray();

                // --- model values               
                double[] yDataModelled = new double[xData.length];

                for (int i = 0; i < yDataModelled.length; i++) {
                    yDataModelled[i] = eval_t2d(xData[i], params);
                }

                // --- jacobian
                double[][] jacobian = new double[yDataModelled.length][params.length];
                for (int j = 0; j < params.length; j++) { // columns

                    double[] params_plus = params.clone();
                    double[] params_minus = params.clone();
                    
                    double eps = 1e-6* (Math.abs(params[j]) + 1.0);
                    params_plus[j]  += eps;
                    params_minus[j] -= eps;

                    for (int i = 0; i < yDataModelled.length; i++) { // rows
                        double eval_plus = eval_t2d(xData[i], params_plus);
                        double eval_minus = eval_t2d(xData[i], params_minus);
                        jacobian[i][j] = (eval_plus-eval_minus)/(2*eps);
                    }
                }

                return new Pair<>(
                    new ArrayRealVector(yDataModelled),
                    new Array2DRowRealMatrix(jacobian)
                );
            }
        };

        /// --- parameter validator
        ParameterValidator validator = new ParameterValidator() {
            public RealVector validate(RealVector params) {
                
                double[] p = params.toArray();

                // cf. eval_t2d()
                p[6] = Math.min(Math.max(p[6], 0), 50); // transition time t1
                p[8] = Math.min(Math.max(p[8], p[6]), 150); // transition time t2
                p[7] = Math.min(Math.max(p[7], 0.5), 15); // should be small
                p[9] = Math.min(Math.max(p[9], 0.5), 30); // should be small

                p[0] = Math.min(Math.max(p[0], 0), 0.5); // int
                p[1] = Math.min(Math.max(p[1], 0), 0.1); // slope : mm / ms

                //double d1 = p[0]+p[1]*p[6];
                //p[2] = Math.min(Math.max(p[2], 0.9*d1), 1.1*d1);
                p[2] = Math.min(Math.max(p[2], 0), 2.5);
                p[3] = Math.min(Math.max(p[3], 0), 0.1); // slope : mm / ms

                //double d2 = p[2]+p[3]*p[8];
                //p[4] = Math.min(Math.max(p[4], 0.9*d2), 1.1*d2); // int
                p[4] = Math.min(Math.max(p[4], 0), 2.5);
                p[5] = Math.min(Math.max(p[5], 0), 0.1); // slope : mm / ms

                // Contrainte de monotonicité : vérifier sur une grille
                double[] testTimes = {0, 10, 20, 50, 100, 150, 200};
                for (int i = 1; i < testTimes.length; i++) {
                    double d_prev = eval_t2d(testTimes[i-1], p);
                    double d_curr = eval_t2d(testTimes[i],   p);
                    if (d_curr < d_prev) {
                        // Ramener vers les valeurs initiales
                        return new ArrayRealVector(initialValues);
                    }
                }
                
                return new ArrayRealVector(p);
            }
        };

        // from ccdb
        //double[] initialValues = {-0.735266, 	0.1015123, 	-0.025992, 	0.0171851, 	2.312859, 	-0.0070255, 	-4.344386, 	5.362382, 	259.71963, 	112.34450};
        //double[] initialValues = {0, 0.02, 1, 0.02, 2, 0.005, 10, 1, 100, 5};

        /// --- least square problem to solve : yDataModelled should be close to yData
        LeastSquaresProblem problem = new LeastSquaresBuilder().
                                        start(initialValues). // initial parameters
                                        model(model). // the model
                                        target(yData). // data to be fit
                                        parameterValidator(validator). // parameter limits
                                        lazyEvaluation(false).
                                        maxEvaluations(5000).
                                        maxIterations(5000).
                                        build();

        LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);

        double[] fittedParams = optimum.getPoint().toArray();

        return fittedParams;

    }

    /**
     * Update the list of wire time2distance parameters in {@link ResultsOverIterations} based on the list of fit paramenters obtained layer by layer
     * @param layer_list list of parameters layer by layer, the first entry is integrated over all layer
     * @param results is a {@link ResultsOverIterations}
     */
    public static void save_time2distance(ArrayList<double[]> layer_list, ResultsOverIterations results) {
        for (int i = 0; i < 576; i++) {
            AhdcWireId identifier = new AhdcWireId(i);
            double[] p = layer_list.get(AhdcWireId.layer2number(identifier.layer));
            results.time2distance.set(i, p);
        }
    }


    /** Number of threads running simultaneously */
    static int nThreads = 20; // 4
    /** Maximum capacity of the queue conataining events */
    static int queue_capacity = 9000; // 150

    /** Frequency of event loggin */
    static int frequency_of_event_login = 100;

    /** Propagator step size */
    static double stepSize = 2; // mm

    /** Used as a mutex for concurrent access */
    static final Object CCDB_LOCK = new Object();

    /** Max number of events to be processed */
    static double nevents_max = Double.MAX_VALUE;

    /**
     * An interface to encode all histogram analyser
     * 
     * <ul>
     * 
     * <li> {@link AhdcAlignmentAnalyser#analyse_global_histograms(Histos, String, int, ResultsOverIterations)} </li>
     * <li> {@link AhdcAlignmentAnalyser#analyse_layer_histograms(Histos, String, int, ResultsOverIterations)} </li>
     * <li> {@link AhdcAlignmentAnalyser#analyse_layer_2D_histograms(Histos, String, int, ResultsOverIterations)} </li>
     * <li> {@link AhdcAlignmentAnalyser#analyse_wire_histograms(Histos, String, int, ResultsOverIterations)} </li>
     * <li> {@link AhdcAlignmentAnalyser#analyse_wire_2D_histograms(Histos, String, int, ResultsOverIterations)} </li>
     * 
     * </u>
     * 
     */
    @FunctionalInterface
    public interface HistoAnalyser {
        void analyse(Histos histos, String outDir, int niter, ResultsOverIterations results);
    }

    /**
     * Angles are in radians. Compute the difference vetween two angles and ensures that the results is between -pi and pi
     * @param angle1
     * @param angle2
     */
    static double diffAngles(double angle1, double angle2) {
        double dphi = Math.atan2(
			Math.sin(angle1 - angle2),
			Math.cos(angle1 - angle2)
		);
        return dphi;
    }


    /**
     * Uncomment the relevant line to run the analysis
     * 
     * @throws IOException 
     * @throws InterruptedException 
     */
    public static void main(String[] args) throws IOException, InterruptedException {

        /// --- Load inputs from options
        fOptions options = new fOptions("-i", "-o", "-ncpu", "-nevts");
        options.LoadOptions(args);
        options.Show();

        /// --- Verify input that files are not empty
        ArrayList<String> inFiles = verify_files(options.GetValues("-i"));
        if (inFiles.size() == 0) {
            System.out.println("Please provide inputs files using the option: -i");
            return;
        }

        /// --- Check that the output dir exists or create a new one
        String outDir = options.GetValue("-o");
        if (outDir.equals("")) {
            System.out.println("\033[1;31m * Please, provide a working directory with the option : -o\033[0m");
            return;
        }
        check_output_dir(outDir);

        /// --- Read ncpu
        String nThreads_str = options.GetValue("-ncpu");
        if (!nThreads_str.equals("")){
            AhdcAlignmentAnalyser.nThreads = Integer.parseInt(nThreads_str);
        }

        /// --- Read nevent_max
        String nevts_str = options.GetValue("-nevts");
        if (!nevts_str.equals("")){
            AhdcAlignmentAnalyser.nevents_max = Integer.parseInt(nevts_str);
        }

        /// --- Read CCDB path
        String ccdb = System.getenv("CCDB_CONNECTION");
        if (ccdb == null) {
            System.out.println("Please define the env variable : CCDB_CONNECTION");
            return;
        } else {
            if (ccdb.equals("mysql://clas12reader@clasdb-farm.jlab.org/clas12") 
                || ccdb.equals("mysql://clas12writer:geom3try@clasdb.jlab.org/clas12")) {
                System.out.println("\033[31;1mPlease set CCDB_CONNECTION to a sqlite file. \033[0m Eg. sqlite:///[/path2file.sqlite]");
                System.out.println("  ...The code will try to write new tables in the CCDB. It is a safety measure to prevent using clas12reader or clas12writer.");
                return;
            }
        }

        /// --- Choose the application
        
        // // 1) Layer alignment
        // layer_alignment(inFiles, outDir, false,
        //     AlertElasticAnalyser::new, 
        //     AhdcAlignmentAnalyser::analyse_global_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_2D_histograms
        // );

        // // 2) Wire alignment
        // wire_alignment(inFiles, outDir, false,
        //     AlertElasticAnalyser::new,
        //     AhdcAlignmentAnalyser::analyse_global_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_2D_histograms,
        //     AhdcAlignmentAnalyser::analyse_wire_histograms
        //     //,AhdcAlignmentAnalyser::analyse_wire_2D_histograms
        // );

        // // 3) time2distance calibration
        // time2distance_calibration(inFiles, outDir, false, 
        //     AlertElasticAnalyser::new,
        //     AhdcAlignmentAnalyser::analyse_global_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_histograms,
        //     AhdcAlignmentAnalyser::analyse_layer_2D_histograms
        // );

        // 4) New ahdc position scan
        new_ahdc_position_scan(inFiles, outDir, true, 
            AlertGoodTrackFinder::new,
            AhdcAlignmentAnalyser::analyse_global_histograms
        );

        //scan_ahdc_position(args, false);
        //new_ahdc_position_scan(inFiles, outDir, true);
        //layer_alignment(inFiles, outDir, true);
        //wire_alignment(inFiles, outDir, true);
        //time2distance_calibration(inFiles, outDir, true);

        // usage
        // time /w/hallb-scshelf2102/clas12/users/touchte/amon/scripts/hipo/run-ahdc-aligner.sh -i /volatile/clas12/touchte/new-translation-table/reconstructed/022712/elastic-filtered/merged/rec_clas_022712.evio.0-834.hipo -o /lustre24/expphy/volatile/clas12/touchte/new-alignment/wire_alignment_following -ncpu 60 > logger.txt 2>&1
    }

}
