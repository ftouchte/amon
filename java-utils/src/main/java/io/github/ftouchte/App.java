package io.github.ftouchte;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.geom.detector.alert.AHDC.AlertDCWireIdentifier;
import org.jlab.geom.prim.Point3D;

import io.github.ftouchte.alignment.AhdcAlignmentAnalyser;
import io.github.ftouchte.alignment.AhdcWireId;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        System.out.println( "Hello World!" );
        test(args);
    }

    public static void test( String[] args ) {
        int Npts = 10;
        double[] vector = new double[Npts];
        Random rand = new Random();
        for (int i = 0; i < vector.length; i++) {
            vector[i] = Math.abs(rand.nextInt()) % 10;
        }

        int[] sortedIndices = AhdcAlignmentAnalyser.sortIndices(vector, false);

        System.out.println("Vector : " + Arrays.toString(vector));
        System.out.println("Sorted indices : " + Arrays.toString(sortedIndices));

        
        
        System.out.println("************* BEFORE");
        AlertDCFactory factory = new AlertDCFactory();
        AlertDCDetector AHDCdet = factory.createDetectorCLAS(new DatabaseConstantProvider());
        printDetector(AHDCdet);

        System.out.println("************* AFTER");
        double[] layer_angles_start = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};
        //double[] layer_angles_start = {180, 180, 180, 180, 180, 180, 180, 180};
        double[] layer_angles_end = {0.881, 1.023, 0.638, 0.905, 0.545, 0.680, 0.499, 0.567};
        factory.setWireCorrectionAngles(AhdcAlignmentAnalyser.layerAngles2WireAngles(layer_angles_start), AhdcAlignmentAnalyser.layerAngles2WireAngles(layer_angles_end));
        AHDCdet = factory.createDetectorCLAS(new DatabaseConstantProvider());
        printDetector(AHDCdet);



    }

    static void printDetector(AlertDCDetector AHDCdet) {
        int Nwires = 10;
        for (int i = 0; i < Nwires; i++) {
            AlertDCWireIdentifier identifier = new AlertDCWireIdentifier(i);
            //int sector = identifier.getSectorId();
            int sector = 1;
            int number = identifier.getLayerId();
            int layer = number % 10;
            int superlayer = number / 10;
            int component = identifier.getComponentId();
            AlertDCWire wire = AHDCdet.getSector(sector).getSuperlayer(superlayer).getLayer(layer).getComponent(component);
            Point3D origin = wire.getLine().origin();
            Point3D end = wire.getLine().end();
            String str = String.format("%2d  %2d  %2d  %2d % 8.4f % 8.4f % 8.4f % 8.4f % 8.4f % 8.4f", sector, superlayer, layer, component, origin.z(), end.z(), origin.x(), origin.y(), end.x(), end.y());
            System.out.println(str);
        }
    } 
}
