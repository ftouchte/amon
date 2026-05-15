package io.github.ftouchte.kalmanFilter;

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.IOException;
import java.io.File;
import java.io.FileOutputStream;
import java.awt.Font;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.ColumnArrangement;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.chart.title.CompositeTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jlab.groot.data.H1F;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import java.util.ArrayList;
import java.util.List;


public class Renderer {


    /**
     * Convert float[] arrat to double[] array
     * @param input
     * @return double[]
     */
    static double[] toDoubleArray(float[] input) {
        double[] data = new double[input.length];
        for (int i = 0; i < input.length; i++) {
            data[i] = input[i];
        }
        return data;
    }

    /**
     * Add a .pdf extension if the filename do not contains it already
     * @param filename
     */
    public static String sanitizeFilenameForPDF(String filename) {
        return sanitizeFilenameForExtension(filename, ".pdf");
    }

    /**
     * Add a .png extension if the filename do not contains it already
     * @param filename
     */
    public static String sanitizeFilenameForPNG(String filename) {
        return sanitizeFilenameForExtension(filename, ".png");
    }

    /**
     * Add 'extension' at the end of the filename if it does not contains it already
     * @param filename
     * @param extension
     */
    public static String sanitizeFilenameForExtension(String filename, String extension) {
        if (filename == null) 
            return null;

        if (filename.toLowerCase().endsWith(extension)) {
            return filename;
        } else {
            return filename + extension;
        }

    }

    /**
     * Create a XYSeries (a sort of XY curve) from the histogram
     * @param h H1F histo
     * @return
     */
    public static XYSeries create_xy_series(H1F h) {
        double[] xValues = h.getXaxis().getBinCenters();
        double[] xLimits = h.getAxis().getLimits();
        double[] count = toDoubleArray(h.getData());

        // Dataset
        XYSeries series = new XYSeries(h.getName());

        for (int i = 0; i < xValues.length; i++) {
            double lower = xLimits[i];
            double upper = xLimits[i+1];
            series.add(lower, count[i]);
            series.add(upper, count[i]);
        }

        return series;
    }

    /**
     * Create a dataset (here XYSeriesCollection) from all histograms
     * This is an utility to plot multiples histograms on the same chart
     * 
     * @param histoList list of H1F histo
     * @return XYSeriesCollection
     */
    public static XYSeriesCollection create_dataset(ArrayList<H1F> histoList) {
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (H1F h : histoList) {
            XYSeries series = create_xy_series(h);
            dataset.addSeries(series);
        }
        return dataset;
    }

    /**
     * Create a dataset (here XYSeriesCollection) from one histo
     * @param h histo
     * @return XYSeriesCollection
     */
    public static XYSeriesCollection create_dataset(H1F h) {
            return new XYSeriesCollection(create_xy_series(h));
    }

    /**
     * Create JFreeChart from dataset (e.g XYSeriesCollection)
     * @param dataset
     * @return
     */
    public static JFreeChart create_chart_from_dataset(XYSeriesCollection dataset, String title, String xtitle, String ytitle) {
        // Create chart
        JFreeChart chart = ChartFactory.createXYLineChart(title, xtitle, ytitle, dataset, PlotOrientation.VERTICAL, true, true, false);

        return chart;
    }

    /**
     * The user can take this piece of code as example to custimize the chart
     * @param chart working chart
     * @return updated chart
     */
    public static JFreeChart customize_chart_by_default(JFreeChart chart) {

        // Supprimer les ombres et bordures du chart
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.WHITE);

        XYPlot plot = (XYPlot) chart.getPlot();

        /// --- Supprimer le fond gris et la grille
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);

        /// --- Bordure noire épaisse autour du canvas (comme ROOT)
        plot.setOutlineVisible(true);
        // plot.setOutlineStroke(new BasicStroke(1.5f));
        // plot.setOutlinePaint(Color.BLACK);

        /// --- Font
        chart.getTitle().setFont(new Font("Times New Roman", Font.BOLD, 14));
        // x axis
        plot.getDomainAxis().setLabelFont(new Font("Times New Roman", Font.PLAIN, 12));
        plot.getDomainAxis().setTickLabelFont(new Font("Times New Roman", Font.PLAIN, 11));
        // y axis
        plot.getRangeAxis().setLabelFont(new Font("Times New Roman", Font.PLAIN, 12));
        plot.getRangeAxis().setTickLabelFont(new Font("Times New Roman", Font.PLAIN, 11));

        /// --- Control the number of ticks
        // axis
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();
        int nTicksX = 10;
        int nTicksY = 5;
        // Min/Max X
        XYSeriesCollection dataset = (XYSeriesCollection) plot.getDataset();
        double xMin = dataset.getDomainLowerBound(false);
        double xMax = dataset.getDomainUpperBound(false);
        // Min/Max Y
        double yMin = dataset.getRangeLowerBound(false);
        double yMax = dataset.getRangeUpperBound(false);
        // Ticks
        xAxis.setTickUnit(new NumberTickUnit((xMax - xMin) / nTicksX));
        yAxis.setTickUnit(new NumberTickUnit((yMax - yMin) / nTicksY));

        /// --- Others
        // Ticks vers l'intérieur comme ROOT
        // xAxis.setTickMarksVisible(true);
        // xAxis.setTickMarkInsideLength(5.0f);   // ticks vers l'intérieur
        // xAxis.setTickMarkOutsideLength(0.0f);  // pas de ticks vers l'extérieur
        // xAxis.setAxisLineVisible(true);
        // xAxis.setAxisLinePaint(Color.BLACK);
        // xAxis.setAxisLineStroke(new BasicStroke(1.5f));

        
        // Ticks vers l'intérieur comme ROOT
        // yAxis.setTickMarksVisible(true);
        // yAxis.setTickMarkInsideLength(5.0f);
        // yAxis.setTickMarkOutsideLength(0.0f);
        // yAxis.setAxisLineVisible(true);
        // yAxis.setAxisLinePaint(Color.BLACK);
        // yAxis.setAxisLineStroke(new BasicStroke(1.5f));

        /// --- Default line style and color rendering
        XYStepRenderer renderer = new XYStepRenderer();
        renderer.setSeriesPaint(0, Color.BLUE); // couleur du contour
        renderer.setSeriesStroke(0, new BasicStroke(1.5f));
        plot.setRenderer(renderer);

        return chart;
    }

    /**
     * Example of code to modify the color or style of the Series (in our case histogram)
     * @param chart working chart
     * @param SeriesId Series identifer
     * @param color
     * @param stroke linestyle
     * @return updated chart
     */
    public static JFreeChart customize_line_rendering(JFreeChart chart, int SeriesId, Color color, BasicStroke stroke) {
        XYPlot plot = (XYPlot) chart.getPlot();
        XYSeriesCollection dataset = (XYSeriesCollection) plot.getDataset();
        XYStepRenderer renderer = (XYStepRenderer) plot.getRenderer();
        int nSeries = dataset.getSeriesCount();
        if (SeriesId < nSeries) {
            if (color != null)
                renderer.setSeriesPaint(SeriesId, color);
            if (stroke != null)
                renderer.setSeriesStroke(SeriesId, stroke);
        }
        return chart;
    }





    public static JFreeChart create_jfreechart(H1F h) {

        // Create chart
        JFreeChart chart = ChartFactory.createXYLineChart(h.getTitle(), h.getTitleX(), h.getTitleY(), create_dataset(h), PlotOrientation.VERTICAL, false, true, false);

        return chart;
    }

    /**
     * Create stat box for an histogram.
     * 
     * Usage: plot.addAnnotation(create_stat_box(h))
     * 
     * @param h histo
     * @return annotation
     */
    public static XYTitleAnnotation create_stat_box(H1F h) {
        // --- Construire la stat box ---
        BlockContainer container = new BlockContainer(new ColumnArrangement());

        TextTitle titleLine   = new TextTitle("Stats", new Font("Monospaced", Font.BOLD, 12));
        TextTitle entriesLine = new TextTitle("Entries : " + h.getEntries(),               new Font("Monospaced", Font.PLAIN, 10));
        TextTitle meanLine    = new TextTitle(String.format("Mean    : %.4f", h.getMean()), new Font("Monospaced", Font.PLAIN, 10));
        TextTitle stdLine     = new TextTitle(String.format("Std     : %.4f", h.getRMS()),new Font("Monospaced", Font.PLAIN, 10));

        titleLine.setPosition(RectangleEdge.TOP);
        entriesLine.setPosition(RectangleEdge.TOP);
        meanLine.setPosition(RectangleEdge.TOP);
        stdLine.setPosition(RectangleEdge.TOP);

        container.add(titleLine);
        container.add(entriesLine);
        container.add(meanLine);
        container.add(stdLine);

        CompositeTitle statBox = new CompositeTitle(container);
        statBox.setBackgroundPaint(new Color(255, 255, 255, 180));
        statBox.setFrame(new BlockBorder(Color.BLACK));

        // --- Positionner dans le canvas ---
        XYTitleAnnotation annotation = new XYTitleAnnotation(
                0.98, 0.98,  // Nord-Est
                statBox,
                RectangleAnchor.TOP_RIGHT
        );

        return annotation;
    }

    static public void save_histogram_as_png(H1F h, int width, int height, String filename) throws IOException {
        JFreeChart chart = create_jfreechart(h);
        save_jfreechart_as_png(chart, filename, width, height);
    }

    static public void save_histogram_as_png(H1F h, int width, int height) throws IOException {
        save_histogram_as_png(h, width, height, sanitizeFilenameForPNG(h.getName()));
    }

    static public void save_jfreechart_as_png(JFreeChart chart, String filename, int width, int height) throws IOException {
        File file = new File(filename);
        ChartUtils.saveChartAsPNG(file, chart, width, height);
    }

    static public void save_jfreechart_as_pdf(JFreeChart chart, String filename, int width, int height) throws IOException, DocumentException {
        Document document = new Document(new Rectangle(width, height));
        PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream(filename));
        document.open();

        PdfContentByte cb = writer.getDirectContent();
        PdfTemplate template = cb.createTemplate(width, height);
        Graphics2D g2 = new PdfGraphics2D(cb, width, height);

        Rectangle2D chartArea = new Rectangle2D.Double(0, 0, width, height);
        chart.draw(g2, chartArea);
        g2.dispose();

        cb.addTemplate(template, 0, 0);
        document.close();
    }

    static public void save_histogram_as_pdf(H1F h, int width, int height, String filename) throws IOException, DocumentException {
        JFreeChart chart = create_jfreechart(h);
        save_jfreechart_as_pdf(chart, filename, width, height);
    }

    static public void save_histogram_as_pdf(H1F h, int width, int height) throws IOException, DocumentException {
        save_histogram_as_pdf(h, width, height, sanitizeFilenameForPDF(h.getName()));
    }



}