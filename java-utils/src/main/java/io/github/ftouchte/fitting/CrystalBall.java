package io.github.ftouchte.fitting;

import org.apache.commons.math3.special.Erf;
import org.jlab.groot.data.GraphErrors;
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

    // public class QueueSide {
    //     static double left = +1.0;
    //     static double right = -1.0;
    // }

    /**
     * The queue of the crystal ball can be located to the left or to the right.
     * 
     * LEFT (standard) : z = (x-µ)/σ <= -α
     * RIGHT :           z = (x-µ)/σ >= α   , α is always positif
     * 
     * To make a unique condition on -α, we can write: RIGHT for -(x-µ)/σ <= -α,
     * so that z = +(x-µ)/σ for LEFT and  z = -(x-µ)/σ for RIGHT
     */
    public enum QueueSide {

        LEFT(+1), // standard
        RIGHT(-1);

        private final int sign;

        QueueSide(int sign) {
            this.sign = sign;
        }

        public int getSign() {
            return sign;
        }
    }

    QueueSide side = QueueSide.LEFT; // per default

    public void setQueueSide(QueueSide _side) {
        this.side = _side;
    }

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
        double z = side.getSign()*(x - mu)/sigma;
        if (z > -alpha) {
            return Math.exp(-0.5*Math.pow(z,2));
        } else { // z <= -alpha, queue side
            return A*Math.pow(B-z, -npower);
        }
    }

    public double evalFit(double x) {
        return amplitude*eval(x);
    }


    /**
     * Not normalised evaluation. Prevent to create new objetc each time. The queue side is set to the LEFT.
     * @param x
     * @param alpha
     * @param npower
     * @param mu
     * @param sigma
     * @return
     */
    public static double eval(double x, double alpha, double npower, double mu, double sigma) {
        // double A = Math.pow(npower/Math.abs(alpha), npower)*Math.exp(-0.5*Math.pow(alpha,2));
        // double B = npower/Math.abs(alpha) - Math.abs(alpha);
        // double z = (x - mu)/sigma;
        // if (z > -alpha) {
        //     return Math.exp(-0.5*Math.pow(z,2));
        // } else {
        //     return A*Math.pow(B-z, -npower);
        // }
        return CrystalBall.eval(x, alpha, npower, mu, sigma, QueueSide.LEFT);
    }

    /**
     * 
     * @param x
     * @param alpha
     * @param npower
     * @param mu
     * @param sigma
     * @param side queue side, see {@link QueueSide}
     * @return
     */
    public static double eval(double x, double alpha, double npower, double mu, double sigma, QueueSide side) {
        double A = Math.pow(npower/Math.abs(alpha), npower)*Math.exp(-0.5*Math.pow(alpha,2));
        double B = npower/Math.abs(alpha) - Math.abs(alpha);
        double z = side.getSign()*(x - mu)/sigma;
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
        System.out.printf ("   side      :  %+d\n", side.getSign());
        System.out.println("   (fit) amplitude  :  " + amplitude);
        System.out.println("   (fit) cost       :  " + cost);
        System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
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

    /** cost of the least squares fit*/
    double cost = -1; 

    /**
     * 
     * @return {@link #amplitude}
     */
    public double getFitAmplitude() {
        return amplitude;
    }

    public void setFitAmplitude(double _amplitude) {
        amplitude = _amplitude;
    }

    /**
     * 
     * @return {@link #cost}
     */
    public double getFitCost() {
        return cost;
    }

    public void setFitCost(double _cost) {
        cost = _cost;
    }




    /**
     * See: https://en.wikipedia.org/wiki/Crystal_Ball_function#/media/File:CrystalBallFunction.svg
     */
    public static void main(String[] args) {

        // basic plot
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
        C3.setQueueSide(CrystalBall.QueueSide.RIGHT);
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
        canvas.save("/lustre24/expphy/volatile/clas12/touchte/alignment/crytal-ball-pdf.pdf");

        C1.print();
        C2.print();
        C3.print();
    }
}
