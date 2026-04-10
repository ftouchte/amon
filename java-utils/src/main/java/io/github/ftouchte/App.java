package io.github.ftouchte;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

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
            vector[i] = Math.abs(rand.nextInt()) % 1000;
        }

        int[] sortedIndices = AhdcAlignmentAnalyser.sortIndices(vector, false);

        System.out.println("Vector : " + Arrays.toString(vector));
        System.out.println("Sorted indices : " + Arrays.toString(sortedIndices));
    }
}
