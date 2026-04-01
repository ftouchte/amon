package io.github.ftouchte;

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

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.service.alert.ALERTEngine;


/**
 * This code is used in parallel of coatjava{branch: ahdc/alignment}
 */
public class AhdcAlignmentAnalyser {

    static class AhdcWireId {
        /** Number between 0 and 575 */
        int num;
        /** Always 1 */
        int sector;
        /** Can be 11, 21, 22, 31, 32, 41, 42, 51 */
        int layer;
        /** component id on a given layer, numeraotation sarting at 1 */
        int component;
        
        /**
         * Ahdc wire id defined with a num ranging from 0 to 576 (excluded)
         * @param _num
         */
        AhdcWireId(int _num) {
            num = _num;
            int[] res = wire2slc(_num);
            sector = res[0];
            layer = res[1];
            component = res[2];
        }

        /**
         * Ahdc wire id defined with sector, layer, component identifiers
         * @param _sector
         * @param _layer
         * @param _component
         */
        AhdcWireId(int _sector, int _layer, int _component) {
            sector = _sector;
            layer = _layer;
            component = _layer;
            num = slc2wire(_sector, _layer, _component);
        }
    }

    /**
     * Use: add entries and complete the constructor and the merge method.
     */
    static class Histos {

        /** 1D residuals LR per layers */
        ArrayList<H1F> h1_residual_LR_per_layers = new ArrayList<>();
        /** 1D residuals per layers */
        ArrayList<H1F> h1_residual_per_layers = new ArrayList<>();
        /** 1D residuals LR per wire */
        ArrayList<H1F> h1_residual_LR_per_wires = new ArrayList<>();

        Histos(int niter) {
            // h1 residual LR per layers
            for (int i = 0; i < 9; i++) {
                H1F h = new H1F("residual-LR-layer-" + number2layer(i) + "itr-" + niter, "residual LR (layer " + number2layer(i) + ")", 100, -3, 3);
                h.setTitleX("layer " + i + ", residual LR (mm)");
                if (i == 0) h.setTitleX("all layers, residual LR (mm)");
                h.setTitleY("count");
                //h.setOptStat(1111);
                h1_residual_LR_per_layers.add(h);
            }
            // h1 residual
            for (int i = 0; i < 9; i++) {
                H1F h = new H1F("residual-layer-" + number2layer(i) + "itr-" + niter, "residual (layer " + number2layer(i) + ")", 100, -3, 3);
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
        }

        void merge(Histos histos) {
            // Per layers
            for (int i = 0; i < 9; i++) {
                // h1 residual LR
                this.h1_residual_LR_per_layers.get(i).add(histos.h1_residual_LR_per_layers.get(i));
                // h1 residual
                this.h1_residual_per_layers.get(i).add(histos.h1_residual_per_layers.get(i));
            }
            // Per wires
            for (int i = 0; i < 576; i++) {
                this.h1_residual_LR_per_wires.get(i).add(histos.h1_residual_LR_per_wires.get(i));
            }

        }
    }

    static class ResultsOverIterations {
        // angles, layer
        double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_angles_sup = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
        double[] layer_angles_inf = {-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0};
        // residuals, layer
        double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_residuals_sup = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_residuals_inf = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        // angles, wire
        double[] wire_angles = new double[576];
        double[] wire_angles_sup = new double[576];
        double[] wire_angles_inf = new double[576];
        //java.util.Arrays.
        // residuals, wire
        double[] wire_residuals = new double[576];
        double[] wire_residuals_sup = new double[576];
        double[] wire_residuals_inf = new double[576];

        // ResultsOverIterations() {
        //     // initialisation
        //     java.util.Arrays.fill(wire_angles, 0.0);
        //     java.util.Arrays.fill(wire_angles_sup, 3.0);
        //     java.util.Arrays.fill(wire_angles_inf, -3.0);
        //     java.util.Arrays.fill(wire_residuals, 0.0);
        //     java.util.Arrays.fill(wire_residuals_sup, 3.0);
        //     java.util.Arrays.fill(wire_residuals_inf, -3.0);
        // }


        // // --- We already have an estimates of the angles inf and sup
        // double[] layer_angles = {-1.49, 3.82, 3.24, -1.54, -1.96, 3.73, 3.30, -1.61}; // attempt 1
        // // double[] layer_angles = {-1.49, 3.82, 3.24, -1.45, -2.01, 3.73, 3.30, -1.61}; // attempt 2
        // double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        // double[] layer_angles_sup = new double[8];
        // double[] layer_residuals_sup = new double[8];
        // double[] layer_angles_inf = new double[8];
        // double[] layer_residuals_inf = new double[8];
        // for (int i = 0; i < 8; i++) {
        //     layer_angles_sup[i] = layer_angles[i] + 0.1*Math.abs(layer_angles[i]);
        //     layer_angles_inf[i] = layer_angles[i] - 0.1*Math.abs(layer_angles[i]);
        //     layer_residuals_sup[i] = 1.0; // positif
        //     layer_residuals_inf[i] = -1.0; // negatif
        // }

        //    layer 11 --> will be rotated by 0.72 deg 
        //     layer 21 --> will be rotated by 1.52 deg 
        //     layer 22 --> will be rotated by 0.97 deg 
        //     layer 31 --> will be rotated by 0.80 deg 
        //     layer 32 --> will be rotated by 0.35 deg 
        //     layer 41 --> will be rotated by 1.36 deg 
        //     layer 42 --> will be rotated by 0.96 deg 
        //     layer 51 --> will be rotated by 0.67 deg 
    }
    
    /** Number of threads running simultaneously */
    static int nThreads = 10; // 4
    /** Maximum capacity of the queue conataining events */
    static int queue_capacity = 100; // 150

    // public static void main(String[] args) {

    //     /// --- Load inputs from options
    //     fOptions options = new fOptions("-i", "-o");
    //     options.LoadOptions(args);
    //     options.Show();

    //     /// --- Verify input files are not empty
    //     ArrayList<String> inFiles = verify_files(options.GetValues("-i"));
    //     if (inFiles.size() == 0) {
    //         System.out.println("Please provide inputs files using the option: -i");
    //         return;
    //     }

    //     /// --- Check that the output dir exists or create a new one
    //     String outDir = options.GetValue("-o");
    //     check_output_dir(outDir);

    //     /// --- Define initial geometry parameters
    //     AlertDCDetector AHDCdet = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

    //     // --- Results over iterations
    //     ResultsOverIterations results = new ResultsOverIterations();

    //     /// --- rotate AHDC detector
    //     System.out.println(" Initial rotation angles : AHDC detector");
    //     doLayerRotations(AHDCdet, results.layer_angles);

    //     /// --- Global observables
    //     GraphErrors g0 = new GraphErrors("sum-squared-resisual-LR-over-iteration");
    //     g0.setTitle("(1/N) #sqrt{#Sum (residual LR)^2}");
    //     g0.setTitleX("iterations");
    //     g0.setTitleY("cost");

    //     /// --- Loop over criteria
    //     int niter = 0;
    //     double value = 1e10;
    //     //while (niter < 12) {
    //     while (value > 2*1e-3) {
    //         niter++;
    //         // run iteration
    //         System.out.println("\033[1;32m ################################ \033[0m");
    //         System.out.println("\033[1;32m # Start iteration : " + niter + "\033[0m");
    //         System.out.println("\033[1;32m ################################ \033[0m");

    //         run(niter, results, inFiles, outDir, AHDCdet, +70);

    //         // running criteria
    //         value = 0;
    //         int N = 0;
    //         for (int i = 0; i < results.layer_residuals.length; i++) {
    //             value += Math.pow(results.layer_residuals[i], 2);
    //             N++;
    //         }
    //         value = Math.sqrt(value) / N;
    //         g0.addPoint(niter, value, 0, 0);
    //         System.out.println("\033[1;31m =======> convergence criteria : " + value + "\033[0m");
            

    //     } // end loop over criteria / nb iterations
    //     EmbeddedCanvas c0 = new EmbeddedCanvas("canvas-sum-squared-resisual-LR-over-iteration-" + niter,"canvas-sum-squared-resisual-LR-over-iteration-" + niter , 1200, 900);
    //     c0.draw(g0);
    //     c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");

    // }

    static void run(int niter, ResultsOverIterations results, ArrayList<String> inFiles, String outDir, AlertDCDetector AHDCdet, double clas_alignment, boolean flag_run_wire_alignment) {
        
        /// --- Parallelizer
        BlockingQueue<DataEvent> queue = new ArrayBlockingQueue<>(queue_capacity);
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);

        DataEvent EVT_POISON = new HipoDataEvent(new Event());

        /// --- Worker
        // Create as many future (thread result) as the number of threads
        //Future<Histos>[] futures = new Future[nThreads];
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
                    if (nevents % 1000 == 0) {
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
        check_output_dir(outDir + "/layers/");
        String layerOutDir = outDir + "/layers/iter-" + niter;
        check_output_dir(layerOutDir);
        EmbeddedCanvas c = new EmbeddedCanvas(1500, 1200);
        c.divide(3, 3);
        GraphErrors g_layer = new GraphErrors();
        for (int i = 0; i < global_histos.h1_residual_LR_per_layers.size(); i++) {
            H1F h = global_histos.h1_residual_LR_per_layers.get(i);
            //double mean = h.getBinContent(h.getMaximumBin());
            double peak = h.getxAxis().getBinCenter(h.getMaximumBin());
            double sigma = h.getRMS();
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


        if (!flag_run_wire_alignment) return;
        /// --- Analyse histograms per wire : (4x3)x48
        //ArrayList<EmbeddedCanvas> canvas = new ArrayList<>(); // 48 canvas
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


        // /// --- Undo AHDC rotation before applying new rotation angles to prevent accumulation
        // undoLayerRotations(AHDCdet, results.layer_angles);

        // /// --- Update rotation angles
        // computeNewLayerAngles(niter, results);

        // /// --- Rotate AHDC detector
        // doLayerRotations(AHDCdet, results.layer_angles);
        

    }

    static void fill_histos(Histos histos, DataEvent event, AlertDCDetector AHDCdet, AlertElasticAnalyser analyser, ALERTEngine alertEngine) {
        if (analyser.hasElasticElectron(event)) {
            // // Retrieve the kinematics of the electron
            // ParticleRow electron = analyser.getElectron();

            // Run engine
            alertEngine.processDataEventProjOnly(event, AHDCdet);

            // Data analysis
            DataBank trackBank = event.getBank("AHDC::kftrack");
            DataBank hitBank = event.getBank("AHDC::hits");

            // Loop over tracks
            for (int i = 0; i < trackBank.rows(); i++) {
                int trackid = trackBank.getInt("trackid", i);
                for (int j = 0; j < hitBank.rows(); j++) {
                    if (trackid == hitBank.getInt("trackid", j)) {
                        int layer = 10*hitBank.getByte("superlayer", j) + hitBank.getByte("layer", j);
                        double residual = hitBank.getDouble("residual", j); // mm
                        double residual_LR = hitBank.getDouble("timeOverThreshold", j); // mm
                        int component = hitBank.getInt("wire", j);
                        AhdcWireId wireId = new AhdcWireId(1, layer, component);

                        // Fill histograms per layers
                        histos.h1_residual_LR_per_layers.get(layer2number(layer)).fill(residual_LR);
                        histos.h1_residual_per_layers.get(layer2number(layer)).fill(residual);
                        histos.h1_residual_LR_per_layers.get(0).fill(residual_LR);
                        histos.h1_residual_per_layers.get(0).fill(residual);
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

    static void save9Histo1D(ArrayList<H1F> histos, String outDir, String name, int nIter) {
        EmbeddedCanvas c = new EmbeddedCanvas(1500, 1200);
        c.divide(3, 3);
        int i = 0;
        for (H1F h : histos) {
            c.cd(i);
            //h.setOptStat(1111);
            c.draw(h);
            
            // draw legend
            //PaveText paveStats = new PaveText("")

            i++;
        }
        c.save(outDir + "/" + name + "-iter-" + nIter + ".pdf");
        //c.save(outDir + "/" + name + "-iter-" + nIter + ".png");
    }

    /**
     * @brief Convert the digit-layer (11,21,...,51) to layer number between 1 and 8
     * 
     * @param digit 
     * @return layer number
     */
    static int layer2number(int digit) {
        if      (digit == 11) {
            return 1;
        } 
        else if (digit == 21) {
            return 2;
        } 
        else if (digit == 22) {
            return 3;
        } 
        else if (digit == 31) {
            return 4;
        } 
        else if (digit == 32) {
            return 5;
        } 
        else if (digit == 41) {
            return 6;
        } 
        else if (digit == 42) {
            return 7;
        } 
        else if (digit == 51) {
            return 8;
        } else {
            return 0; // not a layer, can encode all layers
        }
    }

    /**
     * @brief Convert the digit-layer (11,21,...,51) to layer number between 1 and 8
     * 
     * @param digit 
     * @return layer number
     */
    static int number2layer(int num) {
        if      (num == 1) {
            return 11;
        } 
        else if (num == 2) {
            return 21;
        } 
        else if (num == 3) {
            return 22;
        } 
        else if (num == 4) {
            return 31;
        } 
        else if (num == 5) {
            return 32;
        } 
        else if (num == 6) {
            return 41;
        } 
        else if (num == 7) {
            return 42;
        } 
        else if (num == 8) {
            return 51;
        } else {
            return 0; // not a layer, can encode all layers
        }
    }

    /**
     * @brief Convert (sector, layer, component) to a unique wire id (number betwwen 0 and 575)
     * 
     * @param sector (not used)
     * @param layer 
     * @param component 
     * @return unique wire id
     */
    static int slc2wire(int sector, int layer, int component) {
        if (layer == 11) {
            return component - 1;
        } 
        else if (layer == 21) {
            return 47 + component - 1;
        } 
        else if (layer == 22) {
            return 47 + 56 + component - 1;
        } 
        else if (layer == 31) {
            return 47 + 56 + 56 + component - 1;
        } 
        else if (layer == 32) {
            return 47 + 56 + 56 + 72 + component - 1;
        } 
        else if (layer == 41) {
            return 47 + 56 + 56 + 72 + 72 + component - 1;
        } 
        else if (layer == 42) {
            return 47 + 56 + 56 + 72 + 72 + 87 + component - 1;
        } 
        else if (layer == 51) {
            return 47 + 56 + 56 + 72 + 72 + 87 + 87 + component - 1;
        } else {
            return -1; // not a ahdc wire
        }
    }

    /**
     * @brief Convert wire number (number from 0 to 575) to (sector,layer,component) ids
     * 
     * This is the invert operation of  @link slc2wire(int, int, int) @endlink 
     * 
     * @param wire wire number between 0 and  576 (excluded)
     * @return a triplet (sector, layer, component) in int[]
     */
    static int[] wire2slc(int wire) {
        int sector = -1;
        int layer = -1;
        int component = -1;
        if (wire < 47) {
            layer = 11;
            component = wire + 1;
        }
        else if ((47 <= wire) && (wire < 47 + 56)) {
            layer = 21;
            component = wire - 47 + 1;
        }
        else if ((47 + 56 <= wire) && (wire < 47 + 56 + 56)) {
            layer = 22;
            component = wire - 47 - 56 + 1;
        }
        else if ((47 + 56 + 56 <= wire) && (wire < 47 + 56 + 56 + 72)) {
            layer = 31;
            component = wire - 47 - 56 - 56 + 1;
        }
        else if ((47 + 56 + 56 + 72 <= wire) && (wire < 47 + 56 + 56 + 72 + 72)) {
            layer = 32;
            component = wire - 47 - 56 - 56 - 72 + 1;
        }
        else if ((47 + 56 + 56 + 72 + 72 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87)) {
            layer = 41;
            component = wire - 47 - 56 - 56 - 72 - 72 + 1;
        }
        else if ((47 + 56 + 56 + 72 + 72 + 87 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87 + 87)) {
            layer = 42;
            component = wire - 47 - 56 - 56 - 72 - 72 - 87 + 1;
        }
        else { // ((47 + 56 + 56 + 72 + 72 + 87 + 87 <= wire) && (wire < 47 + 56 + 56 + 72 + 72 + 87 + 87 + 99)) {
            layer = 51;
            component = wire - 47 - 56 - 56 - 72 - 72 - 87 - 87 + 1;
        }
        return new int[] {sector, layer, component};
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
            int layer = number2layer(num+1);
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
            int layer = number2layer(num+1);
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
        while (value > 2*1e-3 && niter < 25) {
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
            computeNewLayerAngles(niter, results);

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



    public static void main(String[] args) {
        //scan_ahdc_position(args);
        //layer_alignment(args);
        wire_alignment(args);
    }


}
