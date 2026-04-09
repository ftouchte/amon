package io.github.ftouchte.alignment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.util.Pair;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.groot.base.PadMargins;
import org.jlab.groot.base.TStyle;
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
import org.jlab.service.alert.ALERTEngine;

import io.github.ftouchte.filtering.AlertElasticAnalyser;
import io.github.ftouchte.fitting.CrystalBall;
import io.github.ftouchte.fitting.CrystallBallFitter;
import io.github.ftouchte.utils.fOptions;


/**
 * This code is used in parallel of coatjava{branch: ahdc/alignment}. We use coatjva to access the 
 * AHDC geometry and the processEvent() method of some classes.
 * 
 * @note This class contains analysis code ready to used for:
 * - layer alignment
 * - wire alignment
 * - position scan
 */
public class AhdcAlignmentAnalyser {
    
    
    /** Number of threads running simultaneously */
    static int nThreads = 40; // 4
    /** Maximum capacity of the queue conataining events */
    static int queue_capacity = 7000; // 150

    /** Frequency of event loggin */
    static int frequency_of_event_loggin = 100;

    /**
     * This the main method. It is used to run an iteration of the alignment procedure.
     * This code uses parallelism to reduce the computing time. The user is invited to
     * manage the parameters : {@link #nThreads} and {@link #queue_capacity}
     * 
     * 
     * @param niter Iteration number
     * @param results {@link ResultsOverIterations} where to store data over iteration
     * @param inFiles HIPO files containing all events to be processed
     * @param outDir where to store plots
     * @param AHDCdet current state of the AHDC geometry (after rotation)
     * @param clas_alignment position of the center of CLAS with respect to the center of ALERT. Useful for {@link #scan_ahdc_position(String[])}
     * @param flag_run_wire_alignment a flag to prevent fitting of wire by wire histograms due to statistics limitation
     */
    static void run(int niter, ResultsOverIterations results, ArrayList<String> inFiles, String outDir, AlertDCDetector AHDCdet, double clas_alignment, boolean flag_run_wire_alignment) {
        
        /// --- Parallelizer
        BlockingQueue<DataEvent> queue = new ArrayBlockingQueue<>(queue_capacity);
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);

        DataEvent EVT_POISON = new HipoDataEvent(new Event());

        /// --- Worker
        // Create as many future (thread result) as the number of threads
        List<Future<Histos>> futures = new ArrayList<>();

        for (int i = 0; i < nThreads; i++) {
            futures.add(pool.submit(() -> {
                System.out.println("\033[1;32m **** Create new thread : " + Thread.currentThread().getName() + " \033[0m");
                Histos local_histos = new Histos(niter);
                AlertElasticAnalyser analyser = new AlertElasticAnalyser();
                ALERTEngine alertEngine = new ALERTEngine();
                alertEngine.init();
                alertEngine.set_clas_alignement(clas_alignment);
                int nevents = 0;

                while (true) {
                    DataEvent event = queue.take();

                    nevents++;
                    if (nevents % frequency_of_event_loggin == 0) {
                        System.out.println("\033[1;32m > " + Thread.currentThread().getName() + " : \033[0m" + nevents + " events");
                    }
                
                    if (event == EVT_POISON) break;

                    fill_histos(local_histos, event, AHDCdet, analyser, alertEngine);

                }

                return local_histos;
            }));
        }

        /// --- Producer
        for (String file : inFiles) {
            // Initialise HIPO reader for the file
            HipoReader reader = new HipoReader();
            reader.open(file);
            /// Read bank definitions
            SchemaFactory factory = reader.getSchemaFactory();


            while (reader.hasNext()) {
                // load event
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

        /// --- Analyse histograms per layer
        analyse_layer_histograms(global_histos, outDir, niter, results);
        analyse_layer_2D_histograms(global_histos, outDir, niter, results);

        /// --- Analyse histograms per wire
        if (flag_run_wire_alignment) 
            analyse_wire_histograms(global_histos, outDir, niter, results);

        
       
    }

    static void analyse_layer_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        check_output_dir(outDir + "/layers/");
        String layerOutDir = outDir + "/layers/iter-" + niter;
        check_output_dir(layerOutDir);
        EmbeddedCanvas c = new EmbeddedCanvas(1500, 1200);
        c.divide(3, 3);
        GraphErrors g_layer = new GraphErrors();
        for (int i = 0; i < global_histos.h1_residual_LR_per_layers.size(); i++) {

            H1F h = global_histos.h1_residual_LR_per_layers.get(i);

            /// --- initialize fit parameters
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
            double cbSide = +1.0; // per defaulf, queue à gauche
            if (leftWidth > rightWidth) { // right asymmetry
                alpha0 = leftWidth/sigma;
            } else {
                alpha0 = rightWidth/sigma;
                cbSide = -1.0;
            }

            // perform fit
            CrystallBallFitter cbFitter = new CrystallBallFitter();
            cbFitter.setAlphaParameter(alpha0, 0.7*Math.abs(alpha0), + 1.3*Math.abs(alpha0));
            cbFitter.setNpowerParameter(5, 3, 100); // 3 is a standard value
            cbFitter.setMuParameter(xpeak, xmin, xmax);
            cbFitter.setSigmaParameter(sigma, 0.6*sigma, 1.1*sigma);
            cbFitter.setAmplitudeParameter(amp, 0.8*amp, 1.2*amp);
            cbFitter.setQueueSide(cbSide);

            
            Pair<CrystalBall, GraphErrors> fitResult = cbFitter.fit(h, xmin, xmax);
            CrystalBall cb = fitResult.getFirst();
            System.out.println("* layer + " + i);
            cb.print();
            GraphErrors gr = fitResult.getSecond();
            double mean = cb.getMu();

            // store residual for the current rotation angle 
            if (i > 0) {                      
                results.layer_residuals[i-1] = mean;
                // Add results in the title for partical reason (the groot's version is too old)
                h.setTitle(String.format("iter : %d, mean : %.5f mm, alpha : %.3f deg", niter, mean, results.layer_angles[i-1]));
            }
            if (i == 0) {
                h.setTitle(String.format("mean : %.5f", mean));
            }

            g_layer.addPoint(i, mean, 0, 0);

            // Draw histograms
            c.cd(i);
            c.draw(h);
            c.draw(gr, "same L");
            
        }
        c.save(layerOutDir + "/residual_LR_iter_" + niter + ".pdf");
        System.out.println("residual_LR_iter_" + niter + ".pdf created");

        EmbeddedCanvas c_layer = new EmbeddedCanvas(1200, 900);
        g_layer.setTitleX("layer");
        g_layer.setTitleY("mean residual LR");
        g_layer.setTitle("Mean residual LR versus layer");
        c_layer.draw(g_layer);
        c_layer.save(layerOutDir + "/summary-residual-LR-layer-iter-" + niter + ".pdf");
        System.out.println(layerOutDir + "/summary-residual-LR-layer-iter-" + niter + ".pdf");

    }

    static void analyse_wire_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        // (4x3)x48 = 576
        ArrayList<EmbeddedCanvas> canvas = new ArrayList<>(); // 48 canvas
        for (int i = 0; i < 48; i++) {
            EmbeddedCanvas c0 = new EmbeddedCanvas(1500, 1600);
            c0.divide(3, 4);
            canvas.add(c0);
        }
        // create a new dir
        check_output_dir(outDir + "/wires/");
        String wireOutDir = outDir + "/wires/iter-" + niter;
        check_output_dir(wireOutDir);
        GraphErrors g_wire = new GraphErrors();
        for (int i = 0; i < global_histos.h1_residual_LR_per_wires.size(); i++) {
            H1F h = global_histos.h1_residual_LR_per_wires.get(i);
            //double mean = h.getBinContent(h.getMaximumBin());
            double peak = h.getxAxis().getBinCenter(h.getMaximumBin());
            double sigma = h.getRMS();
            if (h.getIntegral() > 100) {
                F1D func = new F1D("func" + h.getName(), "[a]*gaus(x, [b], [c]) + [d]", peak-1.8*sigma, peak+1.8*sigma);
                //F1D func = new F1D("func" + h.getName(), "[a]*gaus(x, [b], [c]) + [d]", -1.0, 1.0);
                func.setParameter(0, h.getMax());
                func.setParameter(1, h.getMean());
                func.setParameter(2, h.getRMS());
                func.setParameter(3, h.getMin());
                func.setParLimits(0, 0, Math.max(1.5*h.getMax(), 1.0)); // prevent to have a trivial interval : lim_min = lim_max
                func.setParLimits(1, -3, 3);
                func.setParLimits(2, 0, 3);
                func.setParLimits(3, 0, Math.max(1.5*h.getMin(),1.0)); // prevent to have a trivial interval : lim_min = lim_max
                func.setLineColor(2);
                func.setLineWidth(2);
                h.fit(func);
                // Retrieve fit results
                double mean = func.getParameter(1);
                double width = func.getParameter(2);

                // store residual for the current rotation angle
                results.wire_residuals[i] = mean;
                h.setTitle(String.format("iter : %d, mean : %.5f mm, alpha : %.3f deg", niter, mean, results.wire_angles[i]));
                g_wire.addPoint(i, mean, 0, 0);
            } else {
                double mean = 0; // set at zero to prevent weird fit // double mean = h.getMean();
                results.wire_residuals[i] = mean;
                g_wire.addPoint(i, mean, 0, 0);
                //h.setTitle("Not enough entry, fit not performed");
                h.setTitle(String.format("(No fit) iter : %d, mean : %.5f mm, alpha : %.3f deg", niter, mean, results.wire_angles[i]));
            }
            

            // draw historgrams
            int canvas_num = i / 12; // between 0 and 47
            int canvas_frame = i % 12; // bbetween 0 and 11
            canvas.get(canvas_num).cd(canvas_frame);
            canvas.get(canvas_num).draw(h);
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
    }

    static void analyse_layer_2D_histograms(Histos global_histos, String outDir, int niter, ResultsOverIterations results) {
        // same outdir as for 1D histograms
        check_output_dir(outDir + "/layers/");
        String layerOutDir = outDir + "/layers/iter-" + niter;
        check_output_dir(layerOutDir);

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
                
                /// --- do a gaussian fit for simplicity
                double peak = h.getxAxis().getBinCenter(h.getMaximumBin());
                double sigma = h.getRMS();
                F1D func = new F1D("func" + h.getName(), "[a]*gaus(x, [b], [c]) + [d]", peak-1.4*sigma, peak+1.4*sigma);
                //F1D func = new F1D("func" + h.getName(), "[a]*gaus(x, [b], [c]) + [d]", -1.0, 1.0);
                func.setParameter(0, h.getMax());
                func.setParameter(1, h.getMean());
                func.setParameter(2, h.getRMS());
                func.setParameter(3, h.getMin());
                func.setParLimits(0, 0, Math.max(1.5*h.getMax(), 1.0)); // prevent to have a trivial interval : lim_min = lim_max
                func.setParLimits(1, -3, 3);
                func.setParLimits(2, 0, 3);
                func.setParLimits(3, 0, Math.max(1.5*h.getMin(),1.0)); // prevent to have a trivial interval : lim_min = lim_max
                func.setLineColor(2);
                func.setLineWidth(2);
                h.fit(func);
                // Retrieve fit results
                double mean = func.getParameter(1);
                double width = func.getParameter(2);

                // store results
                if (width < 0.5){
                    gr.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, width);
                    gr_bis.addPoint(h2.getXAxis().getBinCenter(bin), mean, 0, 0);
                }
            } // end loop over bin

            // now fit the graph with a strait line
            WeightedObservedPoints observedPoints = new WeightedObservedPoints();
            for (int pt = 1; pt < gr_bis.getDataSize(0)-1; pt++) { // the 2 end points are excluded
                observedPoints.add(gr_bis.getDataX(pt), gr_bis.getDataY(pt));
            }
            //observedPoints.add(i, niter);
            PolynomialCurveFitter polFitter = PolynomialCurveFitter.create(2);
            double[] params = polFitter.fit(observedPoints.toList());

            F1D func = new F1D("strait-fit", String.format("%f + %f*x", params[0], params[1]), -16, 16);
            func.setLineColor(2);
            func.setLineWidth(4);
            
            // draw
            c.cd(i);
            c.getPad(i).setPalette("kBird");
            c.getPad(i).getAxisFrame().setDrawAxisZ(false);
            h2_initial.setTitle(String.format("slope : %f, constant : %f", params[1], params[0]));
            c.draw(h2_initial);
            //c.draw(gr, "same P");
            c.draw(func, "same L");
            c.draw(gr_bis, "same P");
            c.getPad(i).getAxisY().setRange(-0.5, 0.5);
            c.getPad(i).getAxisFrame().setDrawAxisZ(false); // again

            

        
        } // end loop over h2

        c.save(layerOutDir + "/corr_vz_residual_LR_iter_" + niter + ".pdf");
        System.out.println("corr_vz_residual_LR_iter_" + niter + ".pdf created");
    }


    static void fill_histos(Histos histos, DataEvent event, AlertDCDetector AHDCdet, AlertElasticAnalyser analyser, ALERTEngine alertEngine) {
        if (analyser.hasElasticElectron(event)) {
            // // Retrieve the kinematics of the electron
            // ParticleRow electron = analyser.getElectron();

            // Run engine
            //alertEngine.processDataEventProjOnly(event, AHDCdet);
            alertEngine.processDataEvent(event, AHDCdet);

            // Data analysis
            DataBank trackBank = event.getBank("AHDC::kftrack");
            DataBank hitBank = event.getBank("AHDC::hits");

            // Loop over tracks
            for (int i = 0; i < trackBank.rows(); i++) {
                int trackid = trackBank.getInt("trackid", i);
                double vz = trackBank.getFloat("z", i)*0.1; // convert mm to cm
                for (int j = 0; j < hitBank.rows(); j++) {
                    if (trackid == hitBank.getInt("trackid", j)) {
                        int layer = 10*hitBank.getByte("superlayer", j) + hitBank.getByte("layer", j);
                        double residual = hitBank.getDouble("residual", j); // mm
                        double residual_LR = hitBank.getDouble("timeOverThreshold", j); // mm
                        int component = hitBank.getInt("wire", j);
                        AhdcWireId wireId = new AhdcWireId(1, layer, component);

                        // Fill histograms per layers
                        histos.h1_residual_LR_per_layers.get(AhdcWireId.layer2number(layer)).fill(residual_LR);
                        histos.h1_residual_per_layers.get(AhdcWireId.layer2number(layer)).fill(residual);
                        histos.h2_corr_vz_residual_LR_per_layers.get(AhdcWireId.layer2number(layer)).fill(vz, residual_LR);
                        histos.h1_residual_LR_per_layers.get(0).fill(residual_LR);
                        histos.h1_residual_per_layers.get(0).fill(residual);
                        histos.h2_corr_vz_residual_LR_per_layers.get(0).fill(vz, residual_LR);

                        // Fill histos per wires
                        histos.h1_residual_LR_per_wires.get(wireId.num).fill(residual_LR);
                    }
                }
                // double px = trackBank.getFloat("px", i);
                // double py = trackBank.getFloat("py", i);
                // double pz = trackBank.getFloat("pz", i);
                // ParticleRow ahdc_track = new ParticleRow(px*Units.MeV, py*Units.MeV, pz*Units.MeV);
                // double vz = trackBank.getFloat("z",i); // mm
                // double sum_adc = trackBank.getInt("sum_adc",i);
            } // end loop over tracks
        }
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
        if (outDir.equals("")) {
            System.out.println("Output dir not provided.");
            String defaultDir = "/w/hallb-scshelf2102/clas12/users/touchte/alignment";
            System.out.println("Check the existence of the default dir : " + defaultDir);
            if (Files.exists(Path.of(defaultDir)) && Files.isDirectory(Path.of(defaultDir))) {
                System.out.println("default directory found");
                outDir = defaultDir;
                System.out.println("set output dir to : " + defaultDir);
            }
            else {
                System.out.println("default directory is missing, please provide a output directory using the option -o");
                return;
            }
        } else {
            if (Files.exists(Path.of(outDir)) && Files.isDirectory(Path.of(outDir))) {
                System.out.println("Output directory is set to : " + outDir);
            }
            else {
                if (!Files.exists(Path.of(outDir))) {
                    try {
                        Files.createDirectories(Path.of(outDir));
                        System.out.println("Output directory created and set to : " + outDir);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    System.out.println("A file or a repo with the same name already exists. Please provide a new output directory using the option -o");
                    return;
                }
                
            }
        }
    }

    /**
     * New method to update rotation angles
     * @param results
     */
    static void computeNewLayerAngles(ResultsOverIterations results) {
        double[] angles = results.layer_angles;
        double[] residuals = results.layer_residuals;
        double r_max = 0;
        for (int i = 0; i < residuals.length; i++) {
            if (Math.abs(residuals[i]) > r_max)
                r_max = Math.abs(residuals[i]);
        }
        for (int i = 0; i < angles.length; i++) {
            double alphaRad = residuals[i]/AhdcWireId.layerNum2Radius(i+1);
            //angles[i] = angles[i] - Math.toDegrees(alphaRad)*Math.pow(residuals[i]/r_max, 2.0)*0.5;
            angles[i] = angles[i] - 0.5*Math.toDegrees(alphaRad);
        }
    }

    static void computeNewLayerAngles(int niter, ResultsOverIterations results) {
        double[] angles = results.layer_angles;
        double[] angles_sup = results.layer_angles_sup;
        double[] angles_inf = results.layer_angles_inf;
        double[] residuals = results.layer_residuals;
        double[] residuals_sup = results.layer_residuals_sup;
        double[] residuals_inf = results.layer_residuals_inf;
        for (int i = 0; i < angles.length; i++) {
            if (niter == 1) {
                if (residuals[i] > 0) {
                    residuals_sup[i] = residuals[i];
                    angles_sup[i] = angles[i];
                    angles[i] = angles[i] - 1;
                } 
                else if (residuals[i] < 0) {
                    residuals_inf[i] = residuals[i];
                    angles_inf[i] = angles[i];
                    angles[i] = angles[i] + 1;
                }
            }
            else if (niter == 2) {
                if (Math.abs(residuals_sup[i]) > 1e-5) { // extrapolate with the sup
                    double slope = (angles[i]-angles_sup[i])/(residuals[i]-residuals_sup[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 0.09*Math.abs(alpha);
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 0.11*Math.abs(alpha);
                    residuals_sup[i] = 1; // artificial, just need to know that is it positif
                    
                }
                else if (Math.abs(residuals_inf[i]) > 1e-5) { // extrapolate with the inf
                    double slope = (angles[i]-angles_inf[i])/(residuals[i]-residuals_inf[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 0.09*Math.abs(alpha);
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 0.11*Math.abs(alpha);
                    residuals_sup[i] = 1; // artificial, just need to know that is it positif
                }
            } else { // now proceed by dichotomie
                if (residuals[i]*residuals_sup[i] < 0) { // residaul_sup is assumed to be positif
                    double tmp = angles[i];
                    angles[i] = (angles_inf[i] + angles_sup[i])/2;
                    angles_inf[i] = tmp;
                    residuals_inf[i] = residuals[i];
                }
                else if (residuals[i]*residuals_inf[i] < 0) { // residaul_sup is assumed to be negatif
                    double tmp = angles[i];
                    angles[i] = (angles_inf[i] + angles_sup[i])/2;
                    angles_sup[i] = tmp;
                    residuals_sup[i] = residuals[i];
                }
            }
        }
    }

    /** Same logic as {@link #computeNewLayerAngles(int, ResultsOverIterations)} but for wires */
    static void computeNewWireAngles(int niter, ResultsOverIterations results) {
        double[] angles = results.wire_angles;
        double[] angles_sup = results.wire_angles_sup;
        double[] angles_inf = results.wire_angles_inf;
        double[] residuals = results.wire_residuals;
        double[] residuals_sup = results.wire_residuals_sup;
        double[] residuals_inf = results.wire_residuals_inf;
        for (int i = 0; i < angles.length; i++) {
            if (Math.abs(residuals[i]) < 1e-9) { //is alrealdy because because we cannot perform fit due to low statistic on a given wire
                continue;
            }
            if (niter == 1) {
                if (residuals[i] > 0) {
                    residuals_sup[i] = residuals[i];
                    angles_sup[i] = angles[i];
                    angles[i] = angles[i] - 1;
                } 
                else if (residuals[i] < 0) {
                    residuals_inf[i] = residuals[i];
                    angles_inf[i] = angles[i];
                    angles[i] = angles[i] + 1;
                }
            }
            else if (niter == 2) {
                if (Math.abs(residuals_sup[i]) > 1e-5) { // extrapolate with the sup
                    double slope = (angles[i]-angles_sup[i])/(residuals[i]-residuals_sup[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 0.09*Math.abs(alpha);
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 0.11*Math.abs(alpha);
                    residuals_sup[i] = 1; // artificial, just need to know that is it positif
                    
                }
                else if (Math.abs(residuals_inf[i]) > 1e-5) { // extrapolate with the inf
                    double slope = (angles[i]-angles_inf[i])/(residuals[i]-residuals_inf[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 0.09*Math.abs(alpha);
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 0.11*Math.abs(alpha); // I make it no symmetric for the dichotomie
                    residuals_sup[i] = 1; // artificial, just need to know that is it positif
                }
            } else { // now proceed by dichotomie
                if (residuals[i]*residuals_sup[i] < 0) { // residaul_sup is assumed to be positif
                    double tmp = angles[i];
                    angles[i] = (angles_inf[i] + angles_sup[i])/2;
                    angles_inf[i] = tmp;
                    residuals_inf[i] = residuals[i];
                }
                else if (residuals[i]*residuals_inf[i] < 0) { // residaul_sup is assumed to be negatif
                    double tmp = angles[i];
                    angles[i] = (angles_inf[i] + angles_sup[i])/2;
                    angles_sup[i] = tmp;
                    residuals_sup[i] = residuals[i];
                }
            }
        }
    }

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

    /** Same logic as {@link #undoLayerRotations(AlertDCDetector, double[])} but for wires */
    static void undoWireRotations(AlertDCDetector AHDCdet, double[] wire_angles) {
        for (int num = 0; num < 576; num++) {
            AhdcWireId id = new AhdcWireId(num);
            // rotate AHDCdet
            int sl = id.layer / 10;
            int l = id.layer % 10;
            int w = id.component;
            AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getComponent(w); // numbering starts at 1
            wire.rotateZ(Math.toRadians(-wire_angles[num]));
        }
    }

    /** Same logic as {@link {@link #doLayerRotations(AlertDCDetector, double[])}} but for wires */
    static void doWireRotations(AlertDCDetector AHDCdet, double[] wire_angles) {
        for (int num = 0; num < 576; num++) {
            AhdcWireId id = new AhdcWireId(num);
            // rotate AHDCdet
            int sl = id.layer / 10;
            int l = id.layer % 10;
            int w = id.component;
            AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getComponent(w); // numbering starts at 1
            wire.rotateZ(Math.toRadians(wire_angles[num]));
            System.out.println(String.format("   wire  %d (L%dC%02d) --> will be rotated by %.2f deg", id.num, id.layer, id.component, wire_angles[num]));
        }
    }

    /**
     * Layer alignment analysis
     */
    static void layer_alignment(String[] args) {

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
        String outDir = options.GetValue("-o");
        check_output_dir(outDir);

        /// --- Define initial geometry parameters
        AlertDCDetector AHDCdet = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        /// --- rotate AHDC detector
        System.out.println(" Initial rotation angles : AHDC detector");
        doLayerRotations(AHDCdet, results.layer_angles);

        /// --- Global observables
        GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
        g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
        g0.setTitleX("iterations");
        g0.setTitleY("cost");

        /// --- Loop over criteria
        int niter = 0;
        double value = 1e10;
        //while (niter < 12) {
        while (value > 1*1e-3 && niter < 50) {
            niter++;
            // run iteration
            System.out.println("\033[1;32m ################################ \033[0m");
            System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
            System.out.println("\033[1;32m ################################ \033[0m");

            run(niter, results, inFiles, outDir, AHDCdet, +75, false);

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

            /// --- Undo AHDC rotation before applying new rotation angles to prevent accumulation
            undoLayerRotations(AHDCdet, results.layer_angles);

            /// --- Update rotation angles
            //computeNewLayerAngles(niter, results);
            computeNewLayerAngles(results);

            /// --- Rotate AHDC detector
            doLayerRotations(AHDCdet, results.layer_angles);

        } // end loop over criteria / nb iterations
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
        System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

    }

    /**
     * Routine to scan ahdc position with respect to CLAS
     */
    static void scan_ahdc_position(String[] args) {

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

                run(niter, results, inFiles, outDir, AHDCdet, clas_alignment, false);

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
                computeNewLayerAngles(niter, results);

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
     * Wire alignment analysis
     */
    static void wire_alignment(String[] args) {

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
        String outDir = options.GetValue("-o");
        check_output_dir(outDir);

        /// --- Define initial geometry parameters
        AlertDCDetector AHDCdet = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

        /// --- Start with a well know layer alignment
        doLayerRotations(AHDCdet, new double[] {0.756, 1.517, 0.984, 0.807, 0.382, 1.305, 0.975, 0.679});

        // --- Results over iterations
        ResultsOverIterations results = new ResultsOverIterations();

        /// --- rotate AHDC detector
        // System.out.println(" Initial rotation angles : AHDC detector");
        // doWireRotations(AHDCdet, results.wire_angles);

        /// --- Global observables
        GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
        g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
        g0.setTitleX("iterations");
        g0.setTitleY("cost");

        /// --- Loop over criteria
        int niter = 0;
        double value = 1e10;
        //while (niter < 12) {
        while (value > 2*1e-3 && niter < 25) {
            niter++;
            // run iteration
            System.out.println("\033[1;32m ################################ \033[0m");
            System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
            System.out.println("\033[1;32m ################################ \033[0m");

            run(niter, results, inFiles, outDir, AHDCdet, +75, true);

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

            /// --- Undo AHDC rotation before applying new rotation angles to prevent accumulation
            undoWireRotations(AHDCdet, results.wire_angles);

            /// --- Update rotation angles
            computeNewWireAngles(niter, results);

            /// --- Rotate AHDC detector
            doWireRotations(AHDCdet, results.wire_angles);
            
        } // end loop over criteria / nb iterations
        EmbeddedCanvas c0 = new EmbeddedCanvas(1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");
        System.out.println(outDir + "/summary-cost-estimation-over-iterations.pdf");

    }


    /**
     * Uncomment the relevant line to run the analysis
     * 
     * Code to be run: amon/scripts/hipo/run-ahdc-aligner.sh
     */
    public static void main(String[] args) {
        //scan_ahdc_position(args);
        layer_alignment(args);
        //wire_alignment(args);
    }


}
