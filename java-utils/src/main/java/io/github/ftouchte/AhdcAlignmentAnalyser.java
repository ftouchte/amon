package io.github.ftouchte;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.io.IOException;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.service.alert.ALERTEngine;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.PaveText;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;


/**
 * This code is used in parallel of coatjava{branch: ahdc/alignment}
 */
public class AhdcAlignmentAnalyser {
    public static void main(String[] args) {

        /// Load options
        fOptions options = new fOptions("-i", "-o");
        options.LoadOptions(args);
        options.Show();

        /// Verify input files are not empty
        ArrayList<String> inFiles = options.GetValues("-i");
        if (inFiles.size() == 0) {
            System.out.println("Please provide inputs files using the option: -i");
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

        // /// Verify that we can perform rotation here
        // AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(1).getLayer(1).getComponent(1);
        // System.out.println("> Before rotation");
        // wire.getLine().show();
        // System.out.println("> After rotation");
        // wire.rotateZ(Math.toRadians(180));
        // wire.getLine().show();
        //if ((superlayerId+1) % 2 == 0) wire.rotateZ(Math.toRadians(1));
        // if      (superlayerId+1 == 1 && layerId+1 == 1) {
        //     wire.rotateZ(Math.toRadians(0));
        // } 
        // else if (superlayerId+1 == 2 && layerId+1 == 1) {
        //     wire.rotateZ(Math.toRadians(2.55));
        // }
        // else if (superlayerId+1 == 2 && layerId+1 == 2) {
        //     wire.rotateZ(Math.toRadians(2.11));
        // }
        // else if (superlayerId+1 == 3 && layerId+1 == 1) {
        //     wire.rotateZ(Math.toRadians(0));
        // }
        // else if (superlayerId+1 == 3 && layerId+1 == 2) {
        //     wire.rotateZ(Math.toRadians(-0.32));
        // }
        // else if (superlayerId+1 == 4 && layerId+1 == 1) {
        //     wire.rotateZ(Math.toRadians(2.37));
        // }
        // else if (superlayerId+1 == 4 && layerId+1 == 2) {
        //     wire.rotateZ(Math.toRadians(2.13));
        // }
        // else if (superlayerId+1 == 5 && layerId+1 == 1) {
        //     wire.rotateZ(Math.toRadians(-0.21));
        // }
        
        



        /// Loop over criteria
        int nbIterations = 0;
        while (nbIterations < 1) {
            nbIterations++;
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
                    if (nevents % 10000 == 0) System.out.println("***** Event " + nevents );
                    //if (nevents > 10) break;

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
                F1D func = new F1D("func" + h.getName(), "[a]*gaus(x, [b], [c])", peak-2*sigma, peak+2*sigma);
                func.setParameter(0, h.getMax());
                func.setParameter(1, h.getMean());
                func.setParameter(2, h.getRMS());
                func.setParLimits(0, 0, 1.5*h.getMax());
                func.setParLimits(1, -3, 3);
                func.setParLimits(2, 0, 3);
                func.setLineColor(2);
                func.setLineWidth(2);
                h.fit(func);
                // Retrieve fit results
                double mean = func.getParameter(1);
                double width = func.getParameter(2);

                // update the rotatio angle of this layer
                if (layer_num > 0) {
                    // estimate rotation angle
                    double angle = 0; // deg
                    // rotate AHDCdet
                    int layer = number2layer(layer_num);
                    int sl = layer / 10;
                    int l = layer % 10;
                    for (int i = 0; i < AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getNumComponents(); i++) {
                        AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(sl).getLayer(l).getComponent(i+1); // numbering starts at 1
                        wire.rotateZ(Math.toRadians(angle));
                    }
                }
                layer_num++;
            }

            /// --- Save histos
            save9Histo1D(h1_residual_LR_per_layers, outDir, "residual_LR", nbIterations);
            //save9Histo1D(h1_residual_per_layers, outDir, "residual");
            
            



            /// --- Geometry update

        } // end loop over criteria


    }

    static void save9Histo1D(ArrayList<H1F> histos, String outDir, String name, int nIter) {
        TGCanvas c = new TGCanvas("canvas_" + name, name, 1500, 1200);
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

    



    //void saveHisto(H1F)
}
