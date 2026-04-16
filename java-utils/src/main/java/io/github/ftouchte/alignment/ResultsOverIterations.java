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

    // crystal ball v4 iter 15
    // public double[] layer_angles_start = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};
    // public double[] layer_angles_end = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};

    // crystal ball v12 iter 12
    public double[] layer_angles_start = {0.8998, 1.0507, 0.6151, 0.8525, 0.4530, 0.6205, 0.4035, 0.5831};
    public double[] layer_angles_end   = {0.8774, 0.9017, 0.5051, 0.8746, 0.5643, 0.6194, 0.4645, 0.4064};

    public double[] layer_residuals_start = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_end = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    public double[] layer_residuals_slope    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_constant = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // wires
    public double[] wire_angles = new double[576];
    public double[] wire_residuals = new double[576];
}
