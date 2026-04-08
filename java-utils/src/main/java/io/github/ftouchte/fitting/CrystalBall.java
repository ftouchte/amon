package io.github.ftouchte.fitting;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.special.Erf;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.graphics.EmbeddedCanvas;



/**
 * 
 * See: https://www.jlab.org/primex/weekly_meetings/slides_2009_07_17/dmitry/crystalball.html
 * Also: https://en.wikipedia.org/wiki/Crystal_Ball_function
 */
public class CrystalBall {

    /** should be greater than 0, define the point where, the PDF changes from a power law to a Gaussin law */
    double alpha;

    /** power of the tail, n should be greater than 1 */
    double npower; 

    /** mu parameter of the Gaussian part */
    double mu;

    /** standard deviation of the Gaussian part */
    double sigma;

    double A, B, C, D, N;

    /**
     * Create an instance of the Cristal Ball distribution. The evaluation of this PDF is done using {@link #eval(double)}
     * 
     * @param _alpha should be greater than 0, define the point where, the PDF changes from a power law to a Gaussin law
     * @param _npower power of the tail, n should be greater than 1
     * @param _mu mu of the Gaussian part
     * @param _sigma standard deviation of the Gaussian part
     */
    public CrystalBall(double _alpha, double _npower, double _mu, double _sigma) {
        alpha = _alpha;
        npower = _npower;
        mu = _mu;
        sigma = _sigma;

        A = Math.pow(npower/Math.abs(alpha), npower)*Math.exp(-0.5*Math.pow(alpha,2));
        B = npower/Math.abs(alpha) - Math.abs(alpha);
        C = (npower/Math.abs(alpha))*(1/(npower-1))*Math.exp(-0.5*Math.pow(alpha,2));
        D = Math.sqrt(0.5*Math.PI)*(1+Erf.erf(Math.abs(alpha)/Math.sqrt(2)));
        double inverseN = sigma*(C+D);
        N = 1/inverseN;

    }

    /**
     * Evaluation of the PDF. Do not take into account an external factor as the {@link #amplitude}.
     */
    public double normalisedEval(double x) {
        return N*eval(x);
    }

    /**
     * Evaluation without the normalisation constant
     */
    public double eval(double x) {
        double z = (x - mu)/sigma;
        if (z > -alpha) {
            return Math.exp(-0.5*Math.pow(z,2));
        } else {
            return A*Math.pow(B-z, -npower);
        }
    }

    public double evalFit(double x) {
        return amplitude*eval(x);
    }


    /**
     * Not normalised evaluation. Prevent to create new objetc each time.
     * @param x
     * @param alpha
     * @param npower
     * @param mu
     * @param sigma
     * @return
     */
    public static double eval(double x, double alpha, double npower, double mu, double sigma) {
        double A = Math.pow(npower/Math.abs(alpha), npower)*Math.exp(-0.5*Math.pow(alpha,2));
        double B = npower/Math.abs(alpha) - Math.abs(alpha);
        double z = (x - mu)/sigma;
        if (z > -alpha) {
            return Math.exp(-0.5*Math.pow(z,2));
        } else {
            return A*Math.pow(B-z, -npower);
        }
    }

    /** Print parameters */
    public void print() {
        System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
        System.out.println("Cristal Ball distribution");
        System.out.println("   alpha     :  " + alpha);
        System.out.println("   npower    :  " + npower);
        System.out.println("   mu        :  " + mu);
        System.out.println("   sigma     :  " + sigma);
        System.out.println("   (fit) amplitude  :  " + amplitude);
        System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
    }

    /**
     * Find the {@link CrystalBall} parameters that fit the most the data
     * 
     * @param xData
     * @param yData
     */
    public static CrystalBall fit(double[] xData, double[] yData, double[] initialParams) {

        /**
         * Compute a new vector yDataModelled acoording to the model
         */
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
                    yDataModelled[i] = _amplitude*CrystalBall.eval(xData[i], _alpha, _npower, _mu, _sigma);
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

        //double[] initialParams = {0.0, 0.0, 0.0, 0.0, 0.0};

        // leats square problem to solve : yDataModelled should be close to yData
        LeastSquaresProblem problem = new LeastSquaresBuilder().
                                        start(initialParams).
                                        model(model). // the jacobian is not specified
                                        target(yData).
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
        return cb;

    }

    /**
     * Fit histogram
     * @param h
     * @return
     */
    public static Pair<CrystalBall, GraphErrors> fit(H1F h, double[] initialParams) {
        double[] xData = h.getxAxis().getBinCenters();
        double[] yData = new double[xData.length];
        float[] float_yData = h.getData();
        for (int i = 0; i < xData.length; i++) {
            yData[i] = (double) float_yData[i];
        }
        
        GraphErrors gr = new GraphErrors();
        gr.setLineColor(2);
        gr.setMarkerSize(0);

        CrystalBall cb = CrystalBall.fit(xData, yData, initialParams);

        for (int i = 0; i < xData.length; i++) {
            gr.addPoint(xData[i], cb.evalFit(xData[i]), 0, 0);
        }

        return new Pair<>(cb, gr);
    }

    public static Pair<CrystalBall, GraphErrors> fit(H1F h, double[] initialParams, double xmin, double xmax) {
        // initial data
        double[] xData = h.getxAxis().getBinCenters();
        double[] yData = new double[xData.length];
        float[] float_yData = h.getData();
        for (int i = 0; i < xData.length; i++) {
            yData[i] = (double) float_yData[i];
        }

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
        xData = new double[xValue.size()];
        yData = new double[xValue.size()];

        for (int i = 0; i < xData.length; i++) {
            xData[i] = xValue.get(i);
            yData[i] = yValue.get(i);
        }

        
        GraphErrors gr = new GraphErrors();
        gr.setLineColor(2);
        gr.setMarkerSize(0);

        CrystalBall cb = CrystalBall.fit(xData, yData, initialParams);

        for (int i = 0; i < xData.length; i++) {
            gr.addPoint(xData[i], cb.evalFit(xData[i]), 0, 0);
        }
        
        return new Pair<>(cb, gr);
    }

    public double getMu() {
        return mu;
    }

    public double getSigma() {
        return sigma;
    }

    public double getNpower() {
        return npower;
    }

    public double getAlpha() {
        return alpha;
    }

    /** amplitude parameter from. Use use {@link #evalFit(double)} */
    double amplitude = 1;

    public double getFitAmplitude() {
        return amplitude;
    }

    public void setFitAmplitude(double _amplitude) {
        amplitude = _amplitude;
    }




    /**
     * See: https://en.wikipedia.org/wiki/Crystal_Ball_function#/media/File:CrystalBallFunction.svg
     */
    public static void main(String[] args) {

        //////////////////////
        // ---- basic plot
        //////////////////////
        EmbeddedCanvas canvas = new EmbeddedCanvas(1200, 900);
        GraphErrors gr1 = new GraphErrors();
        GraphErrors gr2 = new GraphErrors();
        GraphErrors gr3 = new GraphErrors();

        gr1.setLineColor(1); // black
        gr2.setLineColor(2); // red
        gr3.setLineColor(3); // green

        gr1.setMarkerSize(0);
        gr2.setMarkerSize(0);
        gr3.setMarkerSize(0);

        int Npts = 500;
        double xmin = -10;
        double xmax = 4;
        CrystalBall C1 = new CrystalBall(10, 2, 0, 1);
        CrystalBall C2 = new CrystalBall(1, 3, 0, 1);
        CrystalBall C3 = new CrystalBall(1, 2, 0, 1);
        for (int i = 0; i < Npts; i++) {
            double x = xmin + i*(xmax-xmin)/(Npts-1);
            double y1 = C1.normalisedEval(x);
            double y2 = C2.normalisedEval(x);
            double y3 = C3.normalisedEval(x);
            gr1.addPoint(x, y1, 0, 0);
            gr2.addPoint(x, y2, 0, 0);
            gr3.addPoint(x, y3, 0, 0);
        }
        canvas.draw(gr1, "L");
        canvas.draw(gr2, "same L");
        canvas.draw(gr3, "same L");
        canvas.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball.pdf");

        //////////////////////////
        // ---- example of fit
        //////////////////////////
        double[] xData = gr1.getVectorX().getArray();
        double[] yData = gr1.getVectorY().getArray();
        Random rand = new Random();
        for (int i = 0; i < yData.length; i++) {
            // add noise
            yData[i] = yData[i] + 0.1*Math.abs(yData[i])*rand.nextGaussian();
        }

        CrystalBall C1_fitted = CrystalBall.fit(xData, yData, new double[] { 12.1, 3, 0, 2, 2});

        // plots
        EmbeddedCanvas canvas2 = new EmbeddedCanvas(1200, 900);

        GraphErrors gr1bis = new GraphErrors();
        GraphErrors gr1bis_fitted = new GraphErrors();

        gr1.setLineColor(1); // black
        gr1bis_fitted.setLineColor(2); // red

        gr1.setMarkerSize(0);
        gr1bis.setMarkerSize(4);
        gr1bis_fitted.setMarkerSize(0);

        for (int i = 0; i < yData.length; i++) {
            gr1bis.addPoint(xData[i], yData[i], 0, 0);
            gr1bis_fitted.addPoint(xData[i], C1_fitted.evalFit(xData[i]), 0, 0);
        }

        canvas2.draw(gr1, "L");
        canvas2.draw(gr1bis, "same P");
        canvas2.draw(gr1bis_fitted, "same L");

        canvas2.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball-fit.pdf");



        // printing
        System.out.println("Original distribution");
        C1.print();;
        System.out.println("Fitted distribution");
        C1_fitted.print();

        //////////////////////////
        // ---- example of fit 2
        //////////////////////////       
        {
            xData = new double[] {
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

            yData = new double[] {
                1,2,2,5,4,2,4,4,5,3,4,4,15,9,8,14,15,15,12,12,26,18,21,24,
                22,35,21,31,36,39,38,57,45,50,68,64,81,89,84,97,95,118,104,113,
                153,128,133,175,167,187,166,164,162,178,149,149,120,128,136,139,
                118,101,95,78,69,73,75,57,58,53,31,36,25,26,23,18,24,15,10,18,
                11,8,9,13,7,6,6,8,5,6,4,5,4,3,3,3,3,3,4,1
                };

            C1_fitted = CrystalBall.fit(xData, yData, new double[] { 1.5, 3, 0, 0.43, 180});
            System.out.println("Fitted distribution 2");
            C1_fitted.print();

            gr1bis = new GraphErrors();
            gr1bis_fitted = new GraphErrors();

            gr1.setLineColor(1); // black
            gr1bis_fitted.setLineColor(2); // red

            gr1bis.setMarkerSize(4);
            gr1bis_fitted.setMarkerSize(0);

            for (int i = 0; i < yData.length; i++) {
                gr1bis.addPoint(xData[i], yData[i], 0, 0);
                gr1bis_fitted.addPoint(xData[i], C1_fitted.evalFit(xData[i]), 0, 0);
            }

            canvas2 = new EmbeddedCanvas(1200, 900);
            canvas2.draw(gr1, "L");
            canvas2.draw(gr1bis, "same P");
            canvas2.draw(gr1bis_fitted, "same L");

            canvas2.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball-fit-example-2.pdf");
        } // end example 2
        
        

        
    }
}
