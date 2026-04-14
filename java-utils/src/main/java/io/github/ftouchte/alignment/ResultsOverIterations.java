package io.github.ftouchte.alignment;

/**
 * The alignment procedure is made over iterations. This class is used
 * to store informations over iterations.
 */
public class ResultsOverIterations {
    
    // layers
    //public double[] layer_angles = {0.991, 1.245, 0.708, 1.086, 0.638, 1.066, 0.667, 0.951}; // layer-alignment-v2 (no fit)
    //public double[] layer_angles = {0.885, 1.032, 0.647, 0.924, 0.572, 0.771, 0.540, 0.628}; // crystal ball fit v3 (iter 50)
    public double[] layer_angles = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567}; // crystal ball fit v4 (iter 15)
    //public double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    public double[] layer_angles_start = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};
    public double[] layer_angles_end = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};

    public double[] layer_residuals_start = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_end = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // wires
    public double[] wire_angles = new double[576];
    public double[] wire_residuals = new double[576];
}
