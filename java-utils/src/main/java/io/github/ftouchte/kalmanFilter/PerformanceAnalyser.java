package io.github.ftouchte.kalmanFilter;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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

import io.github.ftouchte.alignment.AhdcWireId;
import io.github.ftouchte.filtering.AlertElasticAnalyser;
import io.github.ftouchte.filtering.AlertTrackSelector;
import io.github.ftouchte.utils.ParticleRow;
import io.github.ftouchte.utils.Units;
import io.github.ftouchte.utils.fOptions;
import kotlin.Unit;

/**
 * Utility to study the perfomance of the Kalman Filter. Especially the dependency of the computing time
 * with the propagation stepper size and the number of iteration
 * 
 * The parallelization of this code is inspired by the one of io.github.ftouchte.alignment.AhdcAlignmentAnalyser.run()
 */
public class PerformanceAnalyser {
    
    /** Number of threads running simultaneously */
    int nThreads = 10; // 4

    /** Maximum capacity of the queue conataining events */
    int queue_capacity = 9000; // 150

    /** Frequency of event loggin */
    int frequency_of_event_login = 100;

    final double clas_alignment = +75 * Units.mm; // mm

    int KF_Niter = 40;
    double stepper_size = 0.5; // mm

    public static void main(String[] args) {
        /// --- Load inputs from options
        fOptions options = new fOptions("-i", "-o", "-ncpu");
        options.LoadOptions(args);
        options.Show();

        ArrayList<String> inFiles = options.GetValues("-i");
        String nThreads_str = options.GetValue("-ncpu");
        int nThreads = 10;
        if (!nThreads_str.equals("")){
            nThreads = Integer.parseInt(nThreads_str);
        }

        PerformanceAnalyser pfAnalyser = new PerformanceAnalyser();
        pfAnalyser.set_nThreads(nThreads);
        // pfAnalyser.set_KF_Niter(40);
        // pfAnalyser.set_step_size(0.5);

        // Test 1 : loop over Niter
        ArrayList<Histos> histos_array = new ArrayList<>(); 
        for (int i = 1; i < 40; i++) {
            pfAnalyser.set_KF_Niter(i);
            Histos h = pfAnalyser.run(inFiles);
            histos_array.add(h);
            System.out.println("# KF Niter : " + i);
            h.print();
        }
    }


    Histos run(ArrayList<String> inFiles) {

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

                Histos local_histos = new Histos();
                AlertElasticAnalyser analyser = new AlertElasticAnalyser();

                // Config ALERT engine
                ALERTEngine alertEngine = new ALERTEngine();
                alertEngine.init();
                alertEngine.set_clas_alignement(clas_alignment);
                alertEngine.set_KF_Niter(KF_Niter);
                alertEngine.setStepSize(stepper_size);
                

                // Config AHDC engine
                AHDCEngine ahdcEngine = new AHDCEngine();
                ahdcEngine.init();

                int nevents = 0;

                while (true) {
                    DataEvent event = queue.take(); // take one Event from the producer

                    nevents++;
                    if (nevents % frequency_of_event_login == 0) {
                        System.out.println("\033[1;32m > " + Thread.currentThread().getName() + " : \033[0m" + nevents + " events");
                    }
                
                    if (event == EVT_POISON) break;

                    // filter good events
                    if (analyser.hasGoodTrack(event)) {
                        ahdcEngine.processDataEvent(event);
                        alertEngine.processDataEvent(event); // the ADC geometry is not relevant in this study
                        local_histos.h1_computing_time.fill(alertEngine.getComputingTime() / 1_000_000.0); // convert ns to ms
                        // process this event
                        processEvent(event, local_histos, analyser);
                    }

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
        Histos global_histos = new Histos();
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

        /// --- Histogram analysis to be placed here
        analyseHistos(global_histos);

        return global_histos;

    }

    void processEvent(DataEvent event, Histos histos, AlertElasticAnalyser analyser) {

        // Data analysis
        DataBank trackBank = event.getBank("AHDC::kftrack");
  
        // Expected track (computed from the lectron kinematics)
        ParticleRow expected_track = analyser.getExpectedTrack();
        double vz0 = expected_track.vz() / Units.cm;
        // double px0 = expected_track.px(Units.GeV);
        // double py0 = expected_track.py(Units.GeV);
        // double pz0 = expected_track.pz(Units.GeV);
        double p0 = expected_track.p(Units.GeV);
        double theta0 = expected_track.theta(Units.rad);
        double phi0 = expected_track.phi(Units.rad);

        histos.h1_p0.fill(p0 / Units.MeV);
        histos.h1_theta0.fill(theta0 / Units.deg);
        histos.h1_phi0.fill(phi0 / Units.deg);
        histos.h1_vz0.fill(vz0 / Units.cm);


        // Reconstructed track
        ArrayList<Integer> trackRows = analyser.getAhdcKFTrackRows();
        
        for (int i = 0; i < trackRows.size(); i++) {
            
            int track_row = trackRows.get(i);
            double vz = trackBank.getFloat("z" , track_row) * Units.mm;
            double px = trackBank.getFloat("px", track_row) * Units.MeV;
            double py = trackBank.getFloat("py", track_row) * Units.MeV;
            double pz = trackBank.getFloat("pz", track_row) * Units.MeV;
            double p = Math.sqrt(Math.pow(px, 2)+Math.pow(py, 2)+Math.pow(pz, 2));
            double theta = Math.acos(pz/p) * Units.rad;
            double phi   = Math.atan2(py, px) * Units.rad;
            if (phi < 0) phi += 2*Math.PI*Units.rad;

            // Fill histos here
            histos.h1_p.fill(p / Units.MeV);
            histos.h1_theta.fill(theta / Units.deg);
            histos.h1_phi.fill(phi / Units.deg);
            histos.h1_vz.fill(vz / Units.cm);

            histos.h1_delta_p.fill(p / Units.MeV - p0 / Units.MeV);
            histos.h1_delta_theta.fill(theta / Units.deg - theta0 / Units.deg);
            histos.h1_delta_phi.fill(phi / Units.deg - phi0 / Units.deg);
            histos.h1_delta_vz.fill(vz / Units.cm - vz0 / Units.cm);


        } // end loop over good tracks
    }

    void analyseHistos(Histos histos) {

    }

    /**
     * Set the number of iterations of the Kalman filter
     * @param _KF_Niter
     */
    void set_KF_Niter(int _KF_Niter) { KF_Niter = _KF_Niter;}

    /**
     * Set the propagtion step size
     * @param _size
     */
    void set_step_size(double _size) {stepper_size = _size;}

    void set_nThreads(int _nThreads) { nThreads = _nThreads;}
}
