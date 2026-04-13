package io.github.ftouchte.fitting;

import java.util.ArrayList;

import java.util.Random;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.graphics.EmbeddedCanvas;

import io.github.ftouchte.fitting.CrystalBall.QueueSide;

public class CrystallBallFitter {

    private FitParameter alphaParameter;
    private FitParameter npowerParameter;
    private FitParameter muParameter;
    private FitParameter sigmaParameter;
    private FitParameter amplitudeParameter;

    public void setAlphaParameter(double _value, double _min, double _max) {
        alphaParameter = new FitParameter(_value, _min, _max);
        alphaParameter.setName("alpha");
    }

    public void setNpowerParameter(double _value, double _min, double _max) {
        npowerParameter = new FitParameter(_value, _min, _max);
        npowerParameter.setName("npower");
    }

    public void setMuParameter(double _value, double _min, double _max) {
        muParameter = new FitParameter(_value, _min, _max);
        muParameter.setName("mu");
    }

    public void setSigmaParameter(double _value, double _min, double _max) {
        sigmaParameter = new FitParameter(_value, _min, _max);
        sigmaParameter.setName("sigma");
    }

    public void setAmplitudeParameter(double _value, double _min, double _max) {
        amplitudeParameter = new FitParameter(_value, _min, _max);
        amplitudeParameter.setName("amplitude");
    }

    public double[] getParameters() {
        return new double[] {alphaParameter.getValue(), npowerParameter.getValue(), muParameter.getValue(), sigmaParameter.getValue(), amplitudeParameter.getValue()};
    }

    public boolean isNotReady() {
        return (alphaParameter == null || npowerParameter == null || muParameter == null || sigmaParameter == null || amplitudeParameter == null);
    }

    private QueueSide side = CrystalBall.QueueSide.LEFT;
    public void setQueueSide(CrystalBall.QueueSide _side) {
        side = _side;
    }

    public CrystalBall.QueueSide getSide() {
        return side;
    }
    
    /**
     * Find the {@link CrystalBall} parameters that fit the most the data
     * 
     * @param xData
     * @param yData
     */
    public CrystalBall fit(double[] xData, double[] yData) {

        if (isNotReady()) {
            System.out.println("Some fit parameters have not been initialised");
            return null;
        }

        // --- Compute a new vector yDataModelled acoording to the model
        MultivariateJacobianFunction model = new MultivariateJacobianFunction() {
            
            public Pair<RealVector, RealMatrix> value(final RealVector parameters) {

                double[] params = parameters.toArray();

                double _alpha = params[0];
                double _npower = params[1];
                double _mu = params[2];
                double _sigma = params[3];
                double _amplitude = params[4];

                // --- model values
                
                //CrystalBall func = new CrystalBall(_alpha, _npower, _mu, _sigma);
                
                double[] yDataModelled = new double[xData.length];

                for (int i = 0; i < yDataModelled.length; i++) {
                    yDataModelled[i] = _amplitude*CrystalBall.eval(xData[i], _alpha, _npower, _mu, _sigma,  side);
                }

                // --- jacobian
                double[][] jacobian = new double[yDataModelled.length][params.length];
                for (int j = 0; j < params.length; j++) { // columns

                    double[] params_plus = params.clone();
                    double[] params_minus = params.clone();
                    
                    double eps = 1e-6* (Math.abs(params[j]) + 1.0);
                    params_plus[j]  += eps;
                    params_minus[j] -= eps;

                    // CrystalBall func_plus = new CrystalBall(params_plus[0], params_plus[1], params_plus[2], params_plus[3]);
                    // CrystalBall func_minus = new CrystalBall(params_minus[0], params_minus[1], params_minus[2], params_minus[3]);

                    for (int i = 0; i < yDataModelled.length; i++) { // rows
                        //jacobian[i][j] = (params_plus[4]*func_plus.eval(xData[i]) - params_minus[4]*func_minus.eval(xData[i]))/(2*eps);
                        double eval_plus = params_plus[4]*CrystalBall.eval(xData[i], params_plus[0], params_plus[1], params_plus[2], params_plus[3]);
                        double eval_minus = params_minus[4]*CrystalBall.eval(xData[i], params_minus[0], params_minus[1], params_minus[2], params_minus[3]);
                        jacobian[i][j] = (eval_plus-eval_minus)/(2*eps);
                    }
                }


                return new Pair<>(
                    new ArrayRealVector(yDataModelled),
                    new Array2DRowRealMatrix(jacobian)
                );
            }
        };

        ParameterValidator validator = new ParameterValidator() {
            public RealVector validate(RealVector params) {
                
                double[] p = params.toArray();

                // alpha ∈ [minVal, maxVal]
                p[0] = Math.min(Math.max(p[0], alphaParameter.getMinValue()), alphaParameter.getMaxValue());

                // npower ∈ [minVal, maxVal]
                p[1] = Math.min(Math.max(p[1], npowerParameter.getMinValue()), npowerParameter.getMaxValue());

                // mu ∈ [minVal, maxVal]
                p[2] = Math.min(Math.max(p[2], muParameter.getMinValue()), muParameter.getMaxValue());

                // sigma ∈ [minVal, maxVal]
                p[3] = Math.min(Math.max(p[3], sigmaParameter.getMinValue()), sigmaParameter.getMaxValue());

                // amplitude ∈ [minVal, maxVal]
                p[4] = Math.min(Math.max(p[4], amplitudeParameter.getMinValue()), amplitudeParameter.getMaxValue());
                
                return new ArrayRealVector(p);
            }
        };




        // least square problem to solve : yDataModelled should be close to yData
        LeastSquaresProblem problem = new LeastSquaresBuilder().
                                        start(this.getParameters()). // initial parameters
                                        model(model). // the model
                                        target(yData). // data to be fit
                                        parameterValidator(validator). // parameter limits
                                        lazyEvaluation(false).
                                        maxEvaluations(1000).
                                        maxIterations(1000).
                                        build();
        
        LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);

        double[] fittedParams = optimum.getPoint().toArray();

        double alpha = fittedParams[0];
        double npower = fittedParams[1];
        double mu = fittedParams[2];
        double sigma = fittedParams[3];
        double amplitude = fittedParams[4];

        // System.out.println("Cristal Ball fit");
        // System.out.println("   alpha     :  " + alpha);
        // System.out.println("   npower    :  " + npower);
        // System.out.println("   mu      :  " + mu);
        // System.out.println("   sigma     :  " + sigma);
        // System.out.println("   integral  :  " + integral);

        CrystalBall cb = new CrystalBall(alpha, npower, mu, sigma);
        cb.setFitAmplitude(amplitude);
        cb.setFitCost(optimum.getCost());
        return cb;

    }

    public CrystalBall fit(double[] xData, double[] yData, double xmin, double xmax) {
        // limited data
        ArrayList<Double> xValue = new ArrayList<>();
        ArrayList<Double> yValue = new ArrayList<>();
        for (int i = 0; i < xData.length; i++) {
            double x = xData[i];
            if (x >= xmin && x <= xmax) {
                xValue.add(x);
                yValue.add(yData[i]);
            }
        }

        // redefinition
        double[] new_xData = new double[xValue.size()];
        double[] new_yData = new double[xValue.size()];

        for (int i = 0; i < xValue.size(); i++) {
            new_xData[i] = xValue.get(i);
            new_yData[i] = yValue.get(i);
        }

        return this.fit(new_xData, new_yData);
    }

    /**
     * Fit histogram
     * @param h
     * @return
     */
    public Pair<CrystalBall, GraphErrors> fit(H1F h) {
        double[] xData = h.getxAxis().getBinCenters();
        double[] yData = float2double(h.getData());

        CrystalBall cb = this.fit(xData, yData);

        GraphErrors gr = new GraphErrors();
        gr.setLineColor(2);
        gr.setMarkerSize(0);

        for (int i = 0; i < xData.length; i++) {
            gr.addPoint(xData[i], cb.evalFit(xData[i]), 0, 0);
        }

        return new Pair<>(cb, gr);
    }

    public Pair<CrystalBall, GraphErrors> fit(H1F h, double xmin, double xmax) {
        // initialisation : convert float[] to double[]
        double[] xData = h.getxAxis().getBinCenters();
        double[] yData = float2double(h.getData());

        CrystalBall cb = this.fit(xData, yData, xmin, xmax);

        GraphErrors gr = new GraphErrors();
        gr.setLineColor(2);
        gr.setMarkerSize(0);

        for (int i = 0; i < xData.length; i++) {
            if (xData[i] >= xmin && xData[i] <= xmax)
                gr.addPoint(xData[i], cb.evalFit(xData[i]), 0, 0);
        }
        
        return new Pair<>(cb, gr);
    }

    private double[] float2double(float[] vector) {
        double[] res = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            res[i] = (double) vector[i];
        }
        return res;
    }

    public void print() {
        System.out.println("> Crsytall Ball Fitter :");
        alphaParameter.print();
        npowerParameter.print();
        muParameter.print();
        sigmaParameter.print();
        amplitudeParameter.print();
    }

    public static void example2(String[] args) {

        double[] xData = new double[] {
            -1.485000, -1.455000, -1.425000, -1.395000, -1.365000, -1.335000,
            -1.305000, -1.275000, -1.245000, -1.215000, -1.185000, -1.155000,
            -1.125000, -1.095000, -1.065000, -1.035000, -1.005000, -0.975000,
            -0.945000, -0.915000, -0.885000, -0.855000, -0.825000, -0.795000,
            -0.765000, -0.735000, -0.705000, -0.675000, -0.645000, -0.615000,
            -0.585000, -0.555000, -0.525000, -0.495000, -0.465000, -0.435000,
            -0.405000, -0.375000, -0.345000, -0.315000, -0.285000, -0.255000,
            -0.225000, -0.195000, -0.165000, -0.135000, -0.105000, -0.075000,
            -0.045000, -0.015000,  0.015000,  0.045000,  0.075000,  0.105000,
            0.135000,  0.165000,  0.195000,  0.225000,  0.255000,  0.285000,
            0.315000,  0.345000,  0.375000,  0.405000,  0.435000,  0.465000,
            0.495000,  0.525000,  0.555000,  0.585000,  0.615000,  0.645000,
            0.675000,  0.705000,  0.735000,  0.765000,  0.795000,  0.825000,
            0.855000,  0.885000,  0.915000,  0.945000,  0.975000,  1.005000,
            1.035000,  1.065000,  1.095000,  1.125000,  1.155000,  1.185000,
            1.215000,  1.245000,  1.275000,  1.305000,  1.335000,  1.365000,
            1.395000,  1.425000,  1.455000,  1.485000
            };

        double[] yData = {
            1,2,2,5,4,2,4,4,5,3,4,4,15,9,8,14,15,15,12,12,26,18,21,24,
            22,35,21,31,36,39,38,57,45,50,68,64,81,89,84,97,95,118,104,113,
            153,128,133,175,167,187,166,164,162,178,149,149,120,128,136,139,
            118,101,95,78,69,73,75,57,58,53,31,36,25,26,23,18,24,15,10,18,
            11,8,9,13,7,6,6,8,5,6,4,5,4,3,3,3,3,3,4,1
            };
        
        CrystallBallFitter fitter = new CrystallBallFitter();

        fitter.setAlphaParameter(1.2, 1, 2);
        fitter.setNpowerParameter(3, 2, 100);
        fitter.setMuParameter(0, -0.1, 0.1);
        fitter.setSigmaParameter(0.4, 0, 1.5);
        fitter.setAmplitudeParameter(180, 170, 200);
        

        double xmin = -0.5;
        double xmax = 0.5;
        CrystalBall cb = fitter.fit(xData, yData, xmin, xmax);
        System.out.println("* Example 2 - fitted distribution");
        cb.print();

        GraphErrors gr = new GraphErrors();
        GraphErrors gr_fitted = new GraphErrors();

        gr.setLineColor(1); // black
        gr_fitted.setLineColor(2); // red

        gr.setMarkerSize(4);
        gr_fitted.setMarkerSize(0);

        for (int i = 0; i < yData.length; i++) {
            gr.addPoint(xData[i], yData[i], 0, 0);
            if (xData[i] >= xmin && xData[i] < xmax)
                gr_fitted.addPoint(xData[i], cb.evalFit(xData[i]), 0, 0);
        }

        EmbeddedCanvas canvas2 = new EmbeddedCanvas(1200, 900);
        canvas2.draw(gr, "same P");
        canvas2.draw(gr_fitted, "same L");

        canvas2.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball-fitter-example-2.pdf");
    }

    public static void example1(String[] args) {

        GraphErrors gr = new GraphErrors();
        gr.setLineColor(1); // black
        gr.setMarkerSize(4);

        // -- data generation
        CrystalBall cb = new CrystalBall(10, 2, 0, 1);
        int Npts = 100;
        double xmin = -10;
        double xmax = 4;
        for (int i = 0; i < Npts; i++) {
            double x = xmin + i*(xmax-xmin)/(Npts-1);
            double y = cb.normalisedEval(x);
            gr.addPoint(x, y, 0, 0);
        }

        // add noise
        GraphErrors gr_error = new GraphErrors();
        //gr_error.setLineColor(1);
        gr_error.setMarkerSize(4);

        double[] xData = gr.getVectorX().getArray();
        double[] yData = gr.getVectorY().getArray();
        Random rand = new Random();
        for (int i = 0; i < yData.length; i++) {
            // add noise
            yData[i] = yData[i] + 0.05*Math.abs(yData[i])*rand.nextGaussian();
            gr_error.addPoint(xData[i], yData[i], 0, 0);
        }

        // fitting
        CrystallBallFitter fitter = new CrystallBallFitter();

        fitter.setAlphaParameter(11, 9, 13);
        fitter.setNpowerParameter(3, 1.3, 5);
        fitter.setMuParameter(0, -0.1, 0.1);
        fitter.setSigmaParameter(0.4, 0, 1.5);
        fitter.setAmplitudeParameter(0.5, 0.3, 1);

        CrystalBall cb_fitted = fitter.fit(xData, yData);
        System.out.println("* Example 1 - fitted distribution");
        cb_fitted.print();

        // plot the fit
        GraphErrors gr_fitted = new GraphErrors();
        gr_fitted.setLineColor(2); // red
        gr_fitted.setMarkerSize(0);

        for (int i = 0; i < xData.length; i++) {
            gr_fitted.addPoint(xData[i], cb_fitted.evalFit(xData[i]), 0, 0);
        }

        EmbeddedCanvas canvas = new EmbeddedCanvas(1200, 900);
        canvas.draw(gr_error, "same P");
        canvas.draw(gr, "same L");
        canvas.draw(gr_fitted, "same L");

        canvas.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball-fitter-example-1.pdf");
    
    }

    public static void main(String[] args) {
        example1(args);
        example2(args);
    }
}
