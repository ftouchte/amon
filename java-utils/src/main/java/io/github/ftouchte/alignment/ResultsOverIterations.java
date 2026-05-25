package io.github.ftouchte.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * The alignment procedure is made over iterations. This class is used
 * to store informations over iterations.
 */
public class ResultsOverIterations {
    
    // layers
    //public double[] layer_angles = {0.991, 1.245, 0.708, 1.086, 0.638, 1.066, 0.667, 0.951}; // layer-alignment-v2 (no fit)
    //public double[] layer_angles = {0.885, 1.032, 0.647, 0.924, 0.572, 0.771, 0.540, 0.628}; // crystal ball fit v3 (iter 50)
    //public double[] layer_angles = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567}; // crystal ball fit v4 (iter 15)
    public double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // old
    // public double[] layer_angles_start = {0.8609, 1.0181, 0.5654, 0.7998, 0.3913, 0.5151, 0.2749, 0.5057};
    // public double[] layer_angles_end   = {0.8412, 0.8157, 0.4084, 0.7939, 0.4747, 0.5086, 0.3518, 0.2534};

    public double[] layer_angles_start = {0.8756, 0.9655, 0.5530, 0.7385, 0.3926, 0.4797, 0.2527, 0.4673};
    public double[] layer_angles_end   = {0.8147, 0.7606, 0.3733, 0.7680, 0.4432, 0.4814, 0.3274, 0.2193};

    public double[] layer_residuals_start = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_end   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    public double[] layer_residuals_slope    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_constant = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    //-----------
    // wires
    //-----------
    public double[] wire_angles = new double[576];
    public double[] wire_residuals = new double[576];

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
            writer.write("# sector, layer, angle_start, angle_end");
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
            writer.write("# sector, layer, wire, angle");
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
    
}
