package io.github.ftouchte.kalmanFilter;

import java.awt.Color;
import java.io.IOException;

import org.jlab.groot.data.H1F;
import org.jlab.groot.math.Dimension1D;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.CategorySeries;
import org.knowm.xchart.VectorGraphicsEncoder;
import org.knowm.xchart.CategorySeries.CategorySeriesRenderStyle;
import org.knowm.xchart.VectorGraphicsEncoder.VectorGraphicsFormat;
import org.knowm.xchart.style.Styler.LegendPosition;

import java.util.ArrayList;
import java.util.List;

public class Renderer {
    
    /**
     * Generate a {@link org.knowm.xchart.CategoryChart#CategoryChart}. The use is free to overwrite the properties of this chart. Example:
     * 
     * CategorySeries series = chart.getSeriesMap("histogram 1");
     * series.setLineColor(Color.RED);
     * 
     * @param h
     * @param height
     * @param width
     * @return CategoryChart object.
     */
    static CategoryChart generate_xchart(H1F h, int width, int height) {

        // Create chart
        CategoryChart chart = 
            new CategoryChartBuilder()
                .width(width)
                .height(height)
                .title(h.getTitle())
                .xAxisTitle(h.getTitleX())
                .yAxisTitle(h.getTitleY())
                .build();
        
        // Customize chart
        chart.getStyler().setLegendPosition(LegendPosition.InsideNW);
        chart.getStyler().setAvailableSpaceFill(.96);
        chart.getStyler().setPlotGridVerticalLinesVisible(false);
        chart.getStyler().setPlotGridHorizontalLinesVisible(false);
        chart.getStyler().setOverlapped(true);

        // Control ticks
        //Dimension1D dimx = new Dimension1D(h.getXaxis().min(), h.getXaxis().max());
        //List<Double> xticks = dimx.getDimensionTicks(10);

        //chart.setCustomXAxisTickLabelsFormatter(null);

        CategorySeries series = chart.addSeries("histogram 1", h.getXaxis().getBinCenters(), toDoubleArray(h.getData()));

        series.setChartCategorySeriesRenderStyle(CategorySeriesRenderStyle.SteppedBar);

        series.setLineColor(Color.BLUE);
        series.setFillColor(new Color(0,0,0,0)); // the transparent coefficient is zero, so no fill color

        return chart;
    }

    static public void save_histogram_as_pdf(H1F h, int width, int height, String filename) throws IOException {
        CategoryChart chart = generate_xchart(h, width, height);
        save_chart_as_pdf(chart, filename);
    }

    static public void save_histogram_as_pdf(H1F h, int width, int height) throws IOException {
        save_histogram_as_pdf(h, width, height, h.getName());
    }


    /**
     * Convert float[] arrat to double[] array
     * @param input
     * @return
     */
    static double[] toDoubleArray(float[] input) {
        double[] data = new double[input.length];
        for (int i = 0; i < input.length; i++) {
            data[i] = input[i];
        }
        return data;
    }

    static public void save_chart_as_pdf(CategoryChart chart, String title) throws IOException {
        VectorGraphicsEncoder.saveVectorGraphic(chart, sanitizeChartTitle(title), VectorGraphicsFormat.PDF);
    }

    public static String sanitizeChartTitle(String title) {
        if (title == null) 
            return null;

        if (title.toLowerCase().endsWith(".pdf")) {
            return title.substring(0, title.length()-4);
        }

        return title;
        
    }

}
