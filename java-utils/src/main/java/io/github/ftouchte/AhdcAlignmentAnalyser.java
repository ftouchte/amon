package io.github.ftouchte;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
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

    static class PairAlphaResidual {
        public Double alpha;
        public Double residual;
        PairAlphaResidual(Double _alpha, Double _residual) {
            alpha = _alpha;
            residual = _residual;
        }
    }

    static class Histos {

        /** 1D residuals LR per layers */
        ArrayList<H1F> h1_residual_LR_per_layers = new ArrayList<>();
        /** 1D residuals per layers */
        ArrayList<H1F> h1_residual_per_layers = new ArrayList<>();

        Histos(int niter) {
            // h1 residual LR
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
        }

        void merge(Histos histos) {
            // Per layers
            for (int i = 0; i < 9; i++) {
                // h1 residual LR
                this.h1_residual_LR_per_layers.get(i).add(histos.h1_residual_LR_per_layers.get(i));
                // h1 residual
                this.h1_residual_per_layers.get(i).add(histos.h1_residual_per_layers.get(i));
            }
        }
    }

    static class ResultsOverIterations {
        // angles
        double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_angles_sup = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
        double[] layer_angles_inf = {-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0};
        // residuals
        double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_residuals_sup = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_residuals_inf = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
    //     TGCanvas c0 = new TGCanvas("canvas-sum-squared-resisual-LR-over-iteration-" + niter,"canvas-sum-squared-resisual-LR-over-iteration-" + niter , 1200, 900);
    //     c0.draw(g0);
    //     c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");

    // }

    static void run(int niter, ResultsOverIterations results, ArrayList<String> inFiles, String outDir, AlertDCDetector AHDCdet, double clas_alignment) {
        
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

        /// --- Analyse histograms
        TGCanvas c = new TGCanvas("canvas-residual-LR-iter-" + niter,"canvas-residual-LR-iter-" + niter , 1500, 1200);
        c.divide(3, 3);
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
            func.setParLimits(0, 0, 1.5*h.getMax());
            func.setParLimits(1, -3, 3);
            func.setParLimits(2, 0, 3);
            func.setParLimits(3, 0, 1.5*h.getMin());
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

            // Draw histograms
            c.cd(i);
            c.draw(h);
            
        }
        c.save(outDir + "/residual_LR_iter_" + niter + ".pdf");

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

                        // Fill histograms
                        histos.h1_residual_LR_per_layers.get(layer2number(layer)).fill(residual_LR);
                        histos.h1_residual_per_layers.get(layer2number(layer)).fill(residual);
                        histos.h1_residual_LR_per_layers.get(0).fill(residual_LR);
                        histos.h1_residual_per_layers.get(0).fill(residual);
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
        TGCanvas c = new TGCanvas("canvas-" + name + "-" + nIter, name, 1500, 1200);
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

    static void computeNewLayerAngles(int niter, ResultsOverIterations results) {
        double[] angles = results.layer_angles;
        double[] angles_sup = results.layer_angles_sup;
        double[] angles_inf = results.layer_angles_inf;
        double[] residuals = results.layer_residuals;
        double[] residuals_sup = results.layer_residuals_sup;
        double[] residuals_inf = results.layer_residuals_inf;
        for (int i = 0; i < 8; i++) {
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
                    double slope = (residuals[i]-residuals_sup[i])/(angles[i]-angles_sup[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 3;
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 3;
                    residuals_sup[i] = 1; // artificial, just need to know that is it positif
                    
                }
                else if (Math.abs(residuals_inf[i]) > 1e-5) { // extrapolate with the inf
                    double slope = (residuals[i]-residuals_inf[i])/(angles[i]-angles_inf[i]);
                    double alpha = slope*(0-residuals[i]) + angles[i]; // here is alpha for residual = 0
                    angles[i] = alpha;
                    angles_inf[i] = alpha - 2;
                    residuals_inf[i] = -1; // artificial, just need to know that is it negatif
                    angles_sup[i] = alpha + 2.5;
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

            run(niter, results, inFiles, outDir, AHDCdet, +75);

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
        TGCanvas c0 = new TGCanvas("canvas-sum-squared-resisual-LR-over-iteration-" + niter,"canvas-sum-squared-resisual-LR-over-iteration-" + niter , 1200, 900);
        c0.draw(g0);
        c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");

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

                run(niter, results, inFiles, outDir, AHDCdet, clas_alignment);

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
            TGCanvas c0 = new TGCanvas("canvas-sum-squared-resisual-LR-over-iteration-" + niter,"canvas-sum-squared-resisual-LR-over-iteration-" + niter , 1200, 900);
            c0.draw(g0);
            c0.save(outDir + "/summary-cost-estimation-over-iterations.pdf");

        } // end loop over steps      
        TGCanvas c0 = new TGCanvas("summary-mean-angle", "summary-mean-angle", 1200, 900);
        c0.draw(g_mean);
        c0.save(outDir0 + "/summary-mean-angle.pdf");
        TGCanvas c1 = new TGCanvas("summary-angle-deviation", "summary-mean-angle", 1200, 900);
        c1.draw(g_dev);
        c0.save(outDir0 + "/summary-angle-deviation.pdf");
        TGCanvas c2 = new TGCanvas("summary-mean-vs-dev", "summary-mean-angle", 1200, 900);
        c2.draw(g_mean_vs_dev);
        c2.save(outDir0 + "/summary-mean-vs-dev.pdf");

    } // end scan ahdc position



    public static void main(String[] args) {
        //scan_ahdc_position(args);
        layer_alignment(args);
    }


}
