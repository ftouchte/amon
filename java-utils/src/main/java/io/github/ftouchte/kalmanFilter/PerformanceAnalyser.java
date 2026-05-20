package io.github.ftouchte.kalmanFilter;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.nio.file.Files;
import java.nio.file.Path;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
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

import com.itextpdf.text.DocumentException;

import io.github.ftouchte.alignment.AhdcWireId;
import io.github.ftouchte.filtering.AlertElasticAnalyser;
import io.github.ftouchte.filtering.AlertTrackSelector;
import io.github.ftouchte.filtering.AlertElasticAnalyser.Mode;
import io.github.ftouchte.kalmanFilter.Renderer.RendererOutputType;
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

    final double clas_alignment = +75; // mm

    int KF_Niter = 40;
    double stepper_size = 0.5; // mm

    public static void main(String[] args) throws IOException, DocumentException {
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

        String workDir = options.GetValue("-o");
        if (workDir.equals("")) {
            System.out.println("\033[1;31m * Please, provide a working directory with the option : -o\033[0m");
            return;
        }

        /// --- Start anlysis
        PerformanceAnalyser pfAnalyser = new PerformanceAnalyser();
        pfAnalyser.set_nThreads(nThreads);
        // pfAnalyser.set_KF_Niter(40);
        // pfAnalyser.set_step_size(0.5);

        /// --- KF Niter
        // int Niter_min = 1;
        // int Niter_max = 2;
        int num_Niter = 2;
        //int[] NiterValues = new int[num_Niter];
        double[] NiterValues = {1, 2, 5, 10, 40};
        //double[] NiterValues = {1, 2};

        /// --- KF stepSize
        // double stepSize_min = 0.5;
        // double stepSize_max = 4;
        //int numStepMinusOne = 10;

        int num_stepSize = 3;
        //double[] stepSizeValues = new double[num_stepSize];
        double[] stepSizeValues = {0.3, 0.5, 1};
        //double[] stepSizeValues = {0.5, 1};

        /// --- Main routine
        /// --- Loop over Niter
        String singleHistoDir = workDir + "/Standalone";
        check_output_dir(singleHistoDir);
        ArrayList<ArrayList<Histos>> HistosCollections = new ArrayList<>();
        //for (int i = Niter_min; i <= Niter_max; i++) {
        for (int i = 0; i < NiterValues.length; i++) {
            int Niter = (int) NiterValues[i];
            pfAnalyser.set_KF_Niter(Niter);
            String kfDir = singleHistoDir + "/KF_Niter" + Niter;
            check_output_dir(kfDir);
            
            ArrayList<Histos> list = new ArrayList<>();
            /// --- Loop over step size
            for (int j = 0; j < stepSizeValues.length; j++) {
                double size = stepSizeValues[j];
                pfAnalyser.set_step_size(size);

                /// --- Run reconstruction
                Histos h = pfAnalyser.run(inFiles, String.format("_KF_Niter_%d_stepSize_%.4f", Niter, size));
                list.add(h);

                // console out
                System.out.printf("# KF Niter : %d  , step size : %f mm\n", Niter, size);
                h.print();
                
                // create directory
                String sizeDir = String.format("%s/stepSize_%.4f", kfDir, size);
                check_output_dir(sizeDir);
                h.save(sizeDir);
            }
            HistosCollections.add(list);
        }

        ArrayList<String> Observables = new ArrayList<>(Arrays.asList("residual", "delta_p", "delta_phi", "delta_theta", "computing_time", "reconstructed_p", "reconstructed_theta", "reconstructed_phi", "reconstructed_vz"));
        String projectionNiterFixedDir = workDir + "/ProjectionNiterFixedStepSizeVarying";
        check_output_dir(projectionNiterFixedDir);
        /// --- Summary plot : Niter fixed, stepSize varies
        for (int i = 0; i < HistosCollections.size(); i++) {
            ArrayList<Histos> list = extractRow(HistosCollections, i);
            for (String obs : Observables) {
                String title = String.format("Evolution of %s with Niter %d over stepSize", obs, (int) NiterValues[i]);
                // plot comobined histograms
                combine_histos(list, obs, title, projectionNiterFixedDir, stepSizeValues, "stepSize");
            }
        }
    
        /// --- Summary plot : stepSize fixed, Niter vary
        String projectionStepSizeFixed = workDir + "/ProjectionStepSizeFixedNiterVarying";
        check_output_dir(projectionStepSizeFixed);
        for (int j = 0; j < HistosCollections.get(0).size(); j++) {
            ArrayList<Histos> list = extractColumn(HistosCollections, j);
            for (String obs : Observables) {
                String title = String.format("Evolution of %s with stepSize %.3f over Niter", obs, stepSizeValues[j]);
                // plot comobined histograms
                combine_histos(list, obs, title, projectionStepSizeFixed, NiterValues, "Niter");
            }
        }




        /// --- Performance study
        // ArrayList<String> Observables = new ArrayList<>(Arrays.asList("residual", "delta_p", "delta_phi", "delta_theta", "computing_time"));
        // for (String obs : Observables) {
        //     ArrayList<Double> Means = new ArrayList<>();
        //     ArrayList<Double> Widths = new ArrayList<>();
        //     for (Histos histos : HistosOverIterations) {
        //         H1F h = histos.getHistogram1D(obs);
        //         double mean = h.getMean();
        //         double width = h.getRMS();
        //         Means.add(Double.valueOf(mean));
        //         Widths.add(Double.valueOf(width));
        //     }
        // }


    }


    Histos run(ArrayList<String> inFiles, String tag) {

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
                analyser.setFilterMode(Mode.IS_ELASTIC);

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
                        //ahdcEngine.processDataEvent(event); // running the ahdc engine can modified the order of the trach rows, some tracks may desappear
                        alertEngine.processDataEvent(event); // the ADC geometry is not relevant in this study
                        local_histos.h1_computing_time.fill(alertEngine.getComputingTime() / 1_000_000.0); // convert ns to ms
                        //System.out.println("diff time " + alertEngine.getComputingTime() / 1_000_000.0);
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
        Histos global_histos = new Histos(tag);
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
        double vz0 = expected_track.vz(); // cm;
        // double px0 = expected_track.px(Units.GeV);
        // double py0 = expected_track.py(Units.GeV);
        // double pz0 = expected_track.pz(Units.GeV);
        double p0 = expected_track.p(Units.GeV);
        double theta0 = expected_track.theta(Units.rad);
        double phi0 = expected_track.phi(Units.rad);
        if (phi0 < 0) phi0 += 2*Math.PI*Units.rad;

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
            double p = Math.sqrt(px*px + py*py + pz*pz);
            double theta = Math.acos(pz/p) * Units.rad;
            double phi   = Math.atan2(py, px) * Units.rad;
            if (phi < 0) phi += 2*Math.PI*Units.rad;

            // Fill histos here
            histos.h1_p.fill(p / Units.MeV);
            histos.h1_theta.fill(theta / Units.deg);
            histos.h1_phi.fill(phi / Units.deg);
            histos.h1_vz.fill(vz / Units.cm);

            histos.h1_delta_p.fill((p0 - p) / Units.MeV);
            histos.h1_delta_theta.fill((theta0 - theta) / Units.deg);
            histos.h1_delta_phi.fill((phi0 - phi) / Units.deg);
            histos.h1_delta_vz.fill((vz0 - vz) / Units.cm);

            DataBank hitBank = event.getBank("AHDC::hits");
            int trackid = trackBank.getInt("trackid", track_row);
            for (int j = 0; j < hitBank.rows(); j++) {
                if (trackid == hitBank.getInt("trackid", j)) {
                    double residual = hitBank.getDouble("residual", j); // mm
                    histos.h1_residual.fill(residual);
                }
            }


        } // end loop over good tracks
    }

    private static ArrayList<Histos> extractRow(ArrayList<ArrayList<Histos>> collection, int i) {
        return i < collection.size() ? collection.get(i) : null;
    }

    private static ArrayList<Histos> extractColumn(ArrayList<ArrayList<Histos>> collection, int j) {
        ArrayList<Histos> list = new ArrayList<>();
        for (ArrayList<Histos> row : collection) {
            list.add(row.get(j));
        }
        return list;
    }

    public static void combine_histos(ArrayList<Histos> collection, String obs, String title, String dir, double[] paramValues, String barName) throws IOException, DocumentException {
            // gather histogram as same name
            ArrayList<H1F> list = new ArrayList<>();
            XYSeries series_mean = new XYSeries("mean");
            XYSeries series_width = new XYSeries("rms");
            for (int i = 0; i < collection.size(); i++) {
                Histos histos = collection.get(i);
                String name = obs + histos.getTag();
                H1F h = histos.getHistogram1D(name);
                list.add(h);

                // extract mean and rms
                double mean = h.getMean();
                double width = h.getRMS();
                series_mean.add(paramValues[i], mean);
                series_width.add(paramValues[i], width);
            }
            // find the reference histogram
            H1F h_ref = null;
            String reconstructed = "reconstructed";
            String expected = "expected";
            if (obs.startsWith(reconstructed)) {
                String refName = expected + obs.substring(reconstructed.length(), obs.length());
                h_ref = collection.get(0).getHistogram1D(refName);
            }
            
            // output
            String filename = title.replace(" ", "_");
            String histoDir = dir + "/combined_histos";
            check_output_dir(histoDir);
            Renderer.generateJetPalette(collection.size());
            JFreeChart chart = Renderer.create_histogram_evolution_with_parameter(list, h_ref, null, title, barName, paramValues);
            Renderer.save_jfreechart(chart, histoDir + "/" + filename, RendererOutputType.PNG);

            // summary plot : mean
            String meanDir = dir + "/mean_histos";
            check_output_dir(meanDir);
            XYSeriesCollection dataset_mean = new XYSeriesCollection();
            dataset_mean.addSeries(series_mean);
            JFreeChart chart_mean = Renderer.create_chart_for_graph(dataset_mean, title, barName, obs);
            XYPlot plot_mean = (XYPlot)chart_mean.getPlot();
            XYLineAndShapeRenderer renderer_mean = (XYLineAndShapeRenderer) plot_mean.getRenderer();
            renderer_mean.setSeriesPaint(0, Color.BLUE);
            renderer_mean.setSeriesLinesVisible(0, true);
            renderer_mean.setSeriesShapesVisible(0, true);
            Renderer.save_jfreechart(chart_mean, meanDir + "/" + filename, RendererOutputType.PNG);

            // summary plot : width
            String widthDir = dir + "/width_histos";
            check_output_dir(widthDir);
            XYSeriesCollection dataset_width = new XYSeriesCollection();
            dataset_width.addSeries(series_width);
            JFreeChart chart_width = Renderer.create_chart_for_graph(dataset_width, title, barName, obs);
            XYPlot plot_width = (XYPlot)chart_width.getPlot();
            XYLineAndShapeRenderer renderer_with = (XYLineAndShapeRenderer) plot_width.getRenderer();
            renderer_with.setSeriesPaint(0, Color.BLUE);
            renderer_with.setSeriesLinesVisible(0, true);
            renderer_with.setSeriesShapesVisible(0, true);
            Renderer.save_jfreechart(chart_width, widthDir + "/" + filename, RendererOutputType.PNG);

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

    static void check_output_dir(String outDir) {
        if (!Files.exists(Path.of(outDir))) {
            try {
                Files.createDirectories(Path.of(outDir));
                System.out.println("Output directory created and set to : " + outDir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            System.out.println("This directory already exist.");
        }
    }

}
