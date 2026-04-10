package io.github.ftouchte.alignment;

/**
 * The alignment procedure is made over iterations. This class is used
 * to store informations over iterations.
 */
public class ResultsOverIterations {
    
    // layers
    public double[] layer_angles = {0.991, 1.245, 0.708, 1.086, 0.638, 1.066, 0.667, 0.951}; // from the alignment without fit
    //public double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // wires
    public double[] wire_angles = new double[576];
    public double[] wire_residuals = new double[576];
}
