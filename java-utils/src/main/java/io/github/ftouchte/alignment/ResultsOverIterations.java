package io.github.ftouchte.alignment;

/**
 * The alignment procedure is made over iterations. This class is used
 * to store informations over iterations.
 */
public class ResultsOverIterations {
    // angles, layer
    public double[] layer_angles = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_angles_sup = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
    public double[] layer_angles_inf = {-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0};
    // residuals, layer
    public double[] layer_residuals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_sup = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public double[] layer_residuals_inf = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // angles, wire
    public double[] wire_angles = new double[576];
    public double[] wire_angles_sup = new double[576];
    public double[] wire_angles_inf = new double[576];

    // residuals, wire
    public double[] wire_residuals = new double[576];
    public double[] wire_residuals_sup = new double[576];
    public double[] wire_residuals_inf = new double[576];
}
