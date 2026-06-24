package io.github.ftouchte.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;

/**
 * The alignment procedure is made over iterations. This class is used
 * to store informations over iterations.
 */
public class ResultsOverIterations {

    ResultsOverIterations() {
        int run = 22712;
		String variation = "default";
        DatabaseConstantProvider cp = new DatabaseConstantProvider(run, variation);
        cp.loadTable("/geometry/alert/ahdc/layer_alignment");
		cp.loadTable("/geometry/alert/ahdc/wire_alignment");
        cp.loadTable("/calibration/alert/ahdc/time_to_distance_wire");

        // Fill layer angles
        for (int i = 0; i < 8; i++) {
            layer_angles_start[i] = cp.getDouble("/geometry/alert/ahdc/layer_alignment/upstream_rotZ", i);
            layer_angles_end[i] = cp.getDouble("/geometry/alert/ahdc/layer_alignment/downstream_rotZ", i);
            layer_residuals_start[i] = 0;
            layer_residuals_end[i] = 0;
            // others
            layer_residuals_slope[i] = 0;
            layer_residuals_constant[i] = 0;
            // just for information
            layer_angles[i] = 0;
            layer_residuals[i] = 0;
        }
        // Fill wire angles
        for (int i = 0; i < 576; i++) {
            wire_angles[i] = cp.getDouble("/geometry/alert/ahdc/wire_alignment/rotZ", i);
            wire_residuals[i] = 0;
        }

        time2distance = new ArrayList<>();
        for (int i = 0; i < 576; i++) {
            double[] params = {
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p1_int", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p1_slope", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p2_int", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p2_slope", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p3_int", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/p3_slope", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/t1_x0", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/t1_width", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/t2_x0", i),
                cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/t2_width", i)
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/z0", i),
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/z1", i),
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/z2", i),
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/extra1", i),
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/extra2", i),
                // cp.getDouble("/calibration/alert/ahdc/time_to_distance_wire/chi2ndf", i)
            };
            time2distance.add(params);
        }

        cp.disconnect();
    }
    
    // layers
    public double[] layer_angles = new double[8];
    public double[] layer_residuals = new double[8];

    // old
    // public double[] layer_angles_start = {0.8609, 1.0181, 0.5654, 0.7998, 0.3913, 0.5151, 0.2749, 0.5057};
    // public double[] layer_angles_end   = {0.8412, 0.8157, 0.4084, 0.7939, 0.4747, 0.5086, 0.3518, 0.2534};

    public double[] layer_angles_start = new double[8];
    public double[] layer_angles_end   = new double[8];

    public double[] layer_residuals_start = new double[8];
    public double[] layer_residuals_end   = new double[8];

    public double[] layer_residuals_slope    = new double[8];
    public double[] layer_residuals_constant = new double[8];

    //-----------
    // wires
    //-----------
    public double[] wire_angles = new double[576];
    public double[] wire_residuals = new double[576];

    // not really used
    public double[] wire_angles_start = new double[576];
    public double[] wire_residuals_start = new double[576];

    public double[] wire_angles_end = new double[576];
    public double[] wire_residuals_end = new double[576];

    public double[] wire_residuals_slope    = new double[576];
    public double[] wire_residuals_constant = new double[576];

    // --------------
    // time2distance
    // --------------
    public ArrayList<double[]> time2distance = null;

    public enum CCDB_TYPE {
        LAYER, // update the layer angles
        WIRE, // update the wire angles
        T2D, // update the tim2distance
    }

    public void save_layer_angles_to_file(String outDir) throws IOException, InterruptedException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "/layer_angles.txt"));
        writer.write("# computed layer angles for the next iteration");
        writer.newLine();
        writer.write("# sector, layer, upstream_rotZ, downstream_rotZ");
        writer.newLine();
        for (int i = 0; i < 8; i++) {
            int sector = 1;
            int layer = AhdcWireId.number2layer(i+1);
            String line = String.format("%2d   %2d   %f   %f", sector, layer, layer_angles_start[i], layer_angles_end[i]);
            writer.write(line);
            writer.newLine();
        }
        writer.close();
    }

    public void save_wire_angles_to_file(String outDir) throws IOException, InterruptedException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "/wire_angles.txt"));
        writer.write("# computed wire angles for the next iteration");
        writer.newLine();
        writer.write("# sector, layer, wire, rotZ");
        writer.newLine();
        for (int i = 0; i < 576; i++) {
            AhdcWireId identifier = new AhdcWireId(i);
            int sector = 1;
            int layer = identifier.layer;
            int wire = identifier.component;
            String line = String.format("%2d   %2d   %2d   %f", sector, layer, wire, wire_angles[i]);
            writer.write(line);
            writer.newLine();
        }
        writer.close();
    }

    public void save_time2distance_to_file(String outDir) throws IOException, InterruptedException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "/time2distance.txt"));
        writer.write("# parameters: (p1_int, p1_slope, p2_int, p2_slope, p3_int, p3_slope, t1_x0, t1_width, t2_x0, t2_width, z0, z1, z2, extra1, extra2, chi2ndf)");
        writer.newLine();
        writer.newLine();
        for (int i = 0; i < 576; i++) {
            AhdcWireId identifier = new AhdcWireId(i);
            double[] p = time2distance.get(i);
            String line = String.format("%2d   %2d   %2d", identifier.sector, identifier.layer, identifier.component);
            for (int j = 0; j < p.length; j++)
                line += "   " + p[j];
            line += "   0.0 0.0 0.0 0.0 0.0 0.0"; // extra parameters
            writer.write(line);
            writer.newLine();
        }
        writer.close();
    }
    
    public void update_ccdb(String outDir, CCDB_TYPE ccdb) throws IOException, InterruptedException {
        String ccdb_env = System.getenv("CCDB_CONNECTION");
        if (ccdb == CCDB_TYPE.LAYER) {
            // Create file
            save_layer_angles_to_file(outDir);
            // Update the CCDB
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_layer.sh", outDir + "/layer_angles.txt");
            pb.environment().put("CCDB_CONNECTION", ccdb_env);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode + "  <--- layer table updated");
        }
        else if (ccdb == CCDB_TYPE.WIRE) {
            // Create file
            save_wire_angles_to_file(outDir);
            // Update CCDB
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_wire.sh", outDir + "/wire_angles.txt");
            pb.environment().put("CCDB_CONNECTION", ccdb_env);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode + "  <--- wire table updated");
        }
        else if (ccdb == CCDB_TYPE.T2D) {
            // Create file
            save_time2distance_to_file(outDir);
            // Update CCDB
            System.out.println("Try to update the T2D ccdb");
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_layer_t2d.sh", outDir + "/time2distance.txt");
            pb.environment().put("CCDB_CONNECTION", ccdb_env);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode + "  <--- time2distance table updated");
        }
    }

    /**
     * Make a screenshot of the layer_angles, wire_angles and time2distance_wire tables
     * <p> Create the files layer_angles.txt, wire_angles.txt, time2distance.txt in the given directory </p>
     * @param outDir
     * @throws IOException
     * @throws InterruptedException
     */
    public void screenshot(String outDir) throws IOException, InterruptedException {
        save_layer_angles_to_file(outDir);
        save_wire_angles_to_file(outDir);
        save_time2distance_to_file(outDir);
    }
    
}
