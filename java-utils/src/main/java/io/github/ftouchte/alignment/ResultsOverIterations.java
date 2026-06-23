package io.github.ftouchte.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.base.ConstantProvider;

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


    public void save(String outDir) throws IOException {
        // layer alignment
        {
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

        // wire alignment
        {
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
    }

    public enum CCDB_TYPE {
        LAYER, // update the layer angles
        WIRE, // update the wire angles
        T2D, // update the tim2distance
        SAVE // save the current state of the ccdb
    }
    
    public void save(String outDir, CCDB_TYPE ccdb) throws IOException, InterruptedException {
        //save(outDir);
        if (ccdb == CCDB_TYPE.LAYER) {
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_layer.sh", outDir + "/layer_angles.txt");
            pb.environment().put("CCDB_CONNECTION", CCDB_CONNECTION);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode);
            System.out.println("\033[1;32m * CCDB updated...\033[0m");
        }
        else if (ccdb == CCDB_TYPE.WIRE) {
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_wire.sh", outDir + "/wire_angles.txt");
            pb.environment().put("CCDB_CONNECTION", CCDB_CONNECTION);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode);
            System.out.println("\033[1;32m * CCDB updated...\033[0m");
        }
        else if (ccdb == CCDB_TYPE.T2D) {
            System.out.println("Try to update the T2D ccdb");
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/update_ccdb_layer_t2d.sh", outDir + "/time2distance.txt");
            pb.environment().put("CCDB_CONNECTION", CCDB_CONNECTION);
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode);
            System.out.println("\033[1;32m * CCDB updated...\033[0m");
        }
        else if (ccdb == CCDB_TYPE.SAVE) {
            ProcessBuilder pb = new ProcessBuilder("bash", "/w/hallb-scshelf2102/clas12/users/touchte/amon/java-utils/src/main/java/io/github/ftouchte/alignment/save_ccdb.sh", outDir);
            Process process = pb.start();
            int exitCode = process.waitFor();
            System.out.println("Script exited with code: " + exitCode);
            System.out.println("\033[1;32m * CCDB saved...\033[0m");
        }
    }

    public static String CCDB_CONNECTION = "sqlite:////volatile/clas12/touchte/new-alignment/test_new_approach_layer/ccdb_2026-05-24.sqlite";

    public static void execute_cmd(String shell, String cmd) throws InterruptedException, IOException {
        ProcessBuilder pb = new ProcessBuilder(shell, "-c", cmd);
        pb.environment().put("CCDB_CONNECTION", CCDB_CONNECTION);
        pb.inheritIO();
        Process process = pb.start();
        int exitCode = process.waitFor();  // attend la fin de l'exécution
        System.out.println("Exit code: " + exitCode + "    (cmd: " + cmd + ")"); // 0 = succès
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        execute_cmd("bash", "ccdb dump /geometry/alert/ahdc/layer_alignment");
        ResultsOverIterations res = new ResultsOverIterations();
        String outDir = "/w/hallb-scshelf2102/clas12/users/touchte/data/simu";
        res.save(outDir);
    }
    
}
