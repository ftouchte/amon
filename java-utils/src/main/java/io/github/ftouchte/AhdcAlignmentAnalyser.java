package io.github.ftouchte;

import java.util.concurrent.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
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

    public static void main(String[] args) {

        /// Load options
        fOptions options = new fOptions("-i", "-o");
        options.LoadOptions(args);
        options.Show();

        /// Verify input files are not empty
        ArrayList<String> inFiles = options.GetValues("-i");
        if (inFiles.size() == 0) {
            System.out.println("Please provide inputs files using the option: -i");
            return;
        }

        /// Check that the output dir exists or create a new one
        String outDir = options.GetValue("-o");
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
                System.out.println("default directory is missing, please provide a output directory using -o. If it does not exist, it will be created");
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
                    System.out.println("A file or a repo with the same name already exists.");
                    return;
                }
                
            }
        }

        /// Initialiase elastic analyser
        AlertElasticAnalyser elasticAnalyser = new AlertElasticAnalyser();

        /// Initialse ATOF, AHDC, ALERT Engine
        // ATOFEngine atofEngine = new ATOFEngine();
        // atofEngine.init();
        // AHDCEngine ahdcEngine = new AHDCEngine();
        // ahdcEngine.init(ModeTrackFinding.AI_Track_Finding);
        ALERTEngine alertEngine = new ALERTEngine();
        alertEngine.init();

        /// --- Define initial geometry parameters
        AlertDCDetector AHDCdet = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

        // --- No memories
        double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_angles_sup = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
        double[] layer_residuals_sup = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] layer_angles_inf = {-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0};
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

        /// rotate AHDC detector
        System.out.println(" Initial rotation angles : AHDC detector");
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

        /// --- Loop over criteria
        int nbIterations = 0;
        while (nbIterations < 35) {
            nbIterations++;
            System.out.println("\033[1;32m\\#######################################\033[0m");
            System.out.println("\033[1;32m\\# Start iteration : " + nbIterations + "\033[0m");
            System.out.println("\033[1;32m\\#######################################\033[0m");
        //     /// --- Initialize histogram
            //H1F h1_residual_LR = new H1F("residual_LR", "residual LR", 100, -3, 3);
            //H1F h1_residual_LR_before = new H1F("residual_LR_before", "residual LR", 100, -3, 3);

            /// Residual per layers
            ArrayList<H1F> h1_residual_LR_per_layers = new ArrayList<>();
            for (int i = 0; i < 9; i++) {
                H1F h = new H1F("residual-LR-layer-" + number2layer(i) + "itr-" + nbIterations, "residual LR (layer " + number2layer(i) + ")", 100, -3, 3);
                h.setTitleX("residual LR (mm)");
                h.setTitleY("count");
                h1_residual_LR_per_layers.add(h);
            }
            ArrayList<H1F> h1_residual_per_layers = new ArrayList<>();
            for (int i = 0; i < 9; i++) {
                H1F h = new H1F("residual-layer-" + number2layer(i) + "itr-" + nbIterations, "residual (layer " + number2layer(i) + ")", 100, -3, 3);
                h.setTitleX("residual (mm)");
                h.setTitleY("count");
                h1_residual_per_layers.add(h);
            }




            /// Loop over file
            for (String file : inFiles) {
                /// Check that the input file exists
                if (Files.exists(Path.of(file))) {
                    System.out.println("> Open file : " + file);
                } else {
                    System.out.println("* " + file + " : error file not found");
                    return;
                }

                /// Initialise HIPO reader for the file
                HipoReader reader = new HipoReader();
                reader.open(file);
                /// Read bank definitions
                SchemaFactory factory = reader.getSchemaFactory();
		        //List<Schema>   schemaList   = factory.getSchemaList();

                Event inEvent = new Event();
                
                int nevents = 0;
                while (reader.hasNext()) { // just check if we have a next event but do not load it
                    /// --- show progress
                    nevents++;
                    if (nevents % 5000 == 0) System.out.println("***** Event " + nevents );
                    if (nevents > 20000) break;

                    /// load event
                    reader.nextEvent(inEvent);
                    DataEvent event = new HipoDataEvent(inEvent, factory);

                    /// Select an event with an elastic electron
                    if (elasticAnalyser.hasElasticElectron(event)) {
                        /// Retrieve the kinematics of the electron
                        ParticleRow electron = elasticAnalyser.getElectron();

                        // /// Before processing
                        // {
                        //     /// Data Analysis
                        //     DataBank trackBank = event.getBank("AHDC::kftrack");
                        //     DataBank hitBank = event.getBank("AHDC::hits");

                        //     for (int i = 0; i < trackBank.rows(); i++) {
                        //         int trackid = trackBank.getInt("trackid", i);
                        //         for (int j = 0; j < hitBank.rows(); j++) {
                        //             if (trackid == hitBank.getInt("trackid", j)) {
                        //                 double residual_LR = hitBank.getDouble("timeOverThreshold", j); // mm
                        //                 // Fill histogram
                        //                 h1_residual_LR_before.fill(residual_LR);
                        //             }
                        //         }
                        //     }
                        // }

                        /// -----------------------------------------
                        /// Run ALERT Engine on this event (all possible tracks, i.e collection of hits are processed)
                        alertEngine.processDataEventProjOnly(event, AHDCdet);
                        /// -----------------------------------------
                        
                        /// Data Analysis
                        DataBank trackBank = event.getBank("AHDC::kftrack");
                        DataBank hitBank = event.getBank("AHDC::hits");

                        for (int i = 0; i < trackBank.rows(); i++) {
                            int trackid = trackBank.getInt("trackid", i);
                            for (int j = 0; j < hitBank.rows(); j++) {
                                if (trackid == hitBank.getInt("trackid", j)) {
                                    int layer = 10*hitBank.getByte("superlayer", j) + hitBank.getByte("layer", j);
                                    double residual = hitBank.getDouble("residual", j); // mm
                                    double residual_LR = hitBank.getDouble("timeOverThreshold", j); // mm
                                    //System.out.println("residual_LR : " + residual_LR);

                                    // Fill histogram
                                    //System.out.println("layer2number(" + layer + ") = " + layer2number(layer));
                                    h1_residual_LR_per_layers.get(layer2number(layer)).fill(residual_LR);
                                    h1_residual_per_layers.get(layer2number(layer)).fill(residual);
                                    h1_residual_LR_per_layers.get(0).fill(residual_LR);
                                    h1_residual_per_layers.get(0).fill(residual);
                                }
                            }
                            double px = trackBank.getFloat("px", i);
                            double py = trackBank.getFloat("py", i);
                            double pz = trackBank.getFloat("pz", i);
                            ParticleRow ahdc_track = new ParticleRow(px*Units.MeV, py*Units.MeV, pz*Units.MeV);
                            double vz = trackBank.getFloat("z",i); // mm
                            double sum_adc = trackBank.getInt("sum_adc",i);
                        } // end loop over track





                    } // end if elastic event

                    
                } // end loop over events



            } // end loop over files

            /// --- Histogram analysis
            int layer_num = 0;
            for (H1F h : h1_residual_LR_per_layers) {
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
                if (layer_num > 0) {                      
                    layer_residuals[layer_num-1] = mean;
                }
                layer_num++;
            } // end histogram analysis

            /// Undo AHDC rotation before applying new rotation angles to prevent accumulation
            undoLayerRotations(AHDCdet, layer_angles);

            computeNewLayerAngles(nbIterations, layer_angles, layer_residuals, layer_angles_sup, layer_residuals_sup, layer_angles_inf, layer_residuals_inf);
            //computeNewLayerAnglesWithMemory(nbIterations, layer_angles, layer_residuals, layer_angles_sup, layer_residuals_sup, layer_angles_inf, layer_residuals_inf);

            /// rotate AHDC detector
            System.out.println(" Niter "  + nbIterations + " : rotate AHDC detector");
            doLayerRotations(AHDCdet, layer_angles);
            

            /// --- Save histos
            save9Histo1D(h1_residual_LR_per_layers, outDir, "residual_LR", nbIterations);
            //save9Histo1D(h1_residual_per_layers, outDir, "residual");
            
            



            /// --- Geometry update

        } // end loop over criteria / nbIterations


    }

    static void save9Histo1D(ArrayList<H1F> histos, String outDir, String name, int nIter) {
        TGCanvas c = new TGCanvas("canvas-" + name + "-" + nIter, name, 1500, 1200);
        c.divide(3, 3);
        int i = 0;
        for (H1F h : histos) {
            c.cd(i);
            h.setOptStat(1111);
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
     * Here angles_inf and angles_sup are not known, so the iterations 1 and 2 are made to give an estimate
     * @param nIter
     * @param angles
     * @param residuals
     * @param angles_sup
     * @param residuals_sup
     * @param angles_inf
     * @param residuals_inf
     */
    static void computeNewLayerAngles(int nIter, double[] angles, double[] residuals, double[] angles_sup, double[] residuals_sup, double[] angles_inf, double[] residuals_inf) {
        for (int i = 0; i < 8; i++) {
            if (nIter == 1) {
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
            else if (nIter == 2) {
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

    /**
     * Here the angles_sup and angles_inf are well defined
     * @param nIter
     * @param angles
     * @param residuals
     * @param angles_sup
     * @param residuals_sup
     * @param angles_inf
     * @param residuals_inf
     */
    static void computeNewLayerAnglesWithMemory(int nIter, double[] angles, double[] residuals, double[] angles_sup, double[] residuals_sup, double[] angles_inf, double[] residuals_inf) {
        for (int i = 0; i < 8; i++) {
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

    



    //void saveHisto(H1F)
}
