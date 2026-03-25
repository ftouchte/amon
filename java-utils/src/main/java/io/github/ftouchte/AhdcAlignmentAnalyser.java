package io.github.ftouchte;

import java.util.ArrayList;
import java.nio.file.Files;
import java.nio.file.Path;

import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriterSorted;
import org.jlab.rec.ahdc.ModeTrackFinding;
import org.jlab.io.base.DataEvent;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.io.base.DataBank;
import org.jlab.service.ahdc.AHDCEngine;
import org.jlab.service.alert.ALERTEngine;
import org.jlab.service.atof.ATOFEngine;


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

        /// Verify that we can perform rotation here
        AlertDCWire wire = AHDCdet.getSector(1).getSuperlayer(1).getLayer(1).getComponent(1);
        System.out.println("> Before rotation");
        wire.getLine().show();
        System.out.println("> After rotation");
        wire.rotateZ(Math.toRadians(180));
        wire.getLine().show();
        



        // /// Loop over criteria
        // while (true) {
        //     /// --- Initialize histogram



        //     /// Loop over file
        //     for (String file : inFiles) {
        //         /// Check that the input file exists
        //         if (Files.exists(Path.of(file))) {
        //             System.out.println("> Input file : " + file);
        //         } else {
        //             System.out.println("* " + file + " : error file not found");
        //             return;
        //         }

        //         /// Initialise HIPO reader for the file
        //         HipoReader reader = new HipoReader();
        //         reader.open(file);

        //         Event inEvent = new Event();
                
        //         int nevents = 0;
        //         while (reader.hasNext()) { // just check if we have a next event but do not load it
        //             /// --- show progress
        //             nevents++;

        //             /// load event
        //             reader.nextEvent(inEvent);
        //             DataEvent event = new HipoDataEvent(inEvent);

        //             /// Select an event with an elastic electron
        //             if (elasticAnalyser.hasElasticElectron(event)) {
        //                 /// Retrieve the kinematics of the electron
        //                 ParticleRow electron = elasticAnalyser.getElectron();
                        
        //                 /// Run ALERT Engine on this event (all possible tracks, i.e collection of hits are processed)
        //                 alertEngine.processDataEventProjOnly(event, AHDCdet);

        //                 /// Data Analysis
        //                 DataBank trackBank = event.getBank("AHDC::kftrack");
        //                 DataBank hitBank = event.getBank("AHDC::hits");

        //                 for (int i = 0; i < trackBank.rows(); i++) {
        //                     int trackid = trackBank.getInt("trackid", i);
        //                     for (int j = 0; j < hitBank.rows(); j++) {
        //                         if (trackid == hitBank.getInt("trackid", j)) {
        //                             double residual = hitBank.getDouble("residual", j); // mm
        //                             double residual_LR = hitBank.getDouble("timeOverThreshold", j); // mm
        //                         }
        //                     }
        //                     double px = trackBank.getFloat("px", i);
        //                     double py = trackBank.getFloat("py", i);
        //                     double pz = trackBank.getFloat("pz", i);
        //                     ParticleRow ahdc_track = new ParticleRow(px*Units.MeV, py*Units.MeV, pz*Units.MeV);
        //                     double vz = trackBank.getFloat("z",i); // mm
        //                     double sum_adc = trackBank.getFloat("sum_adc",i);
        //                 }
        //             }

                    
        //         } // end loop over events



        //     } // end loop over files

        //     /// --- Histogram analysis



        //     /// --- Geometry update

        // } // end loop over criteria


    }
}
