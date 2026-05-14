package io.github.ftouchte;

import java.io.IOException;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;

public class XChartTest {
    public static void main(String[] args) throws IOException {
        double[] xData = new double[]{0.0, 1.0, 2.0};
        double[] yData = new double[]{2.0, 1.0, 0.0};

        
        // Create Chart
        XYChart chart = QuickChart.getChart("Sample Chart", "X", "Y", "y(x)", xData, yData);

        // Save it
        BitmapEncoder.saveBitmap(chart, "./Sample_Chart",BitmapFormat.PNG);

        // Show it
        //new SwingWrapper(chart).displayChart();

    }
}
