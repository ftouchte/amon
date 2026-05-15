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
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.ColumnArrangement;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.chart.title.CompositeTitle;
import org.jfree.chart.title.LegendTitle;
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
        plot.setDomainGridlinesVisible(true);
        plot.setRangeGridlinesVisible(true);

        /// --- Bordure noire épaisse autour du canvas (comme ROOT)
        plot.setOutlineVisible(true);
        plot.setOutlineStroke(new BasicStroke(1.5f));

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
     * @param SeriesId Series identifier
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

        TextTitle titleLine   = new TextTitle("Stats", new Font("Monospaced", Font.BOLD, (int) (1.1*basicFontSize())));
        TextTitle entriesLine = new TextTitle("Entries : " + h.getEntries(),               new Font("Monospaced", Font.PLAIN, (int) basicFontSize()));
        TextTitle meanLine    = new TextTitle(String.format("Mean    : %.4f", h.getMean()), new Font("Monospaced", Font.PLAIN, (int) basicFontSize()));
        TextTitle stdLine     = new TextTitle(String.format("Std     : %.4f", h.getRMS()),new Font("Monospaced", Font.PLAIN, (int) basicFontSize()));

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

    public enum RendererOutputType {
        PDF,
        PNG
    }

    /** Current width of the chart */
    public static int width = 1500;
    /** Current height of the chart */
    public static int height = 1200;
    /**
     * Set dimensions
     * @param _width
     * @param _height
     */
    public static void setDimentions(int _width, int _height) {
        width = _width;
        height = _height;
    }

    /**
     * Access basi font size scaled with the current dimensions of the charts
     * @return
     */
    public static float basicFontSize() {
        return (float) Math.sqrt(width * height) / 45.0f;
    }

    /**
     * Print out a chart with 1 or more histogram
     * @param histos
     * @param filename
     * @param outType
     * @throws IOException
     * @throws DocumentException
     */
    public static void save_overlayed_histogram(ArrayList<H1F> histos, String filename, RendererOutputType outType) throws IOException, DocumentException {
        if (histos.isEmpty() || histos == null) return;
        XYSeriesCollection dataset = create_dataset(histos);
        JFreeChart chart = create_chart_from_dataset(dataset, histos.get(0).getTitle(), histos.get(0).getTitleX(), histos.get(0).getTitleY());
        chart = customize_chart_by_default(chart);
        chart = apply_scalable_fontsize(chart);
        for (int i = 0; i < histos.size(); i++) {
            customize_line_rendering(chart, i, null, null);
        }
        if (histos.size() == 1) {
            XYPlot plot = (XYPlot) chart.getPlot();
            plot.addAnnotation(create_stat_box(histos.get(0)));
            chart.getLegend().setVisible(false);
        }
        if (outType == RendererOutputType.PNG) {
            save_jfreechart_as_png(chart, sanitizeFilenameForPNG(filename));
        }
        else if (outType == RendererOutputType.PDF) {
            save_jfreechart_as_pdf(chart, sanitizeFilenameForPDF(filename));
        }
    }

    /**
     * Print out a chart with 1 histogram
     * @param histo
     * @param filename
     * @param outType
     * @throws IOException
     * @throws DocumentException
     */
    public static void save_histogram(H1F histo, String filename, RendererOutputType outType) throws IOException, DocumentException {
        ArrayList<H1F> list = new ArrayList<>();
        list.add(histo);
        save_overlayed_histogram(list, filename, outType);
    }

    /**
     * Save chart as PNG
     * @param chart
     * @param filename
     * @throws IOException
     */
    static public void save_jfreechart_as_png(JFreeChart chart, String filename) throws IOException {
        File file = new File(filename);
        ChartUtils.saveChartAsPNG(file, chart, width, height);
    }

    /**
     * Save chart as PDF
     * @param chart
     * @param filename
     * @throws IOException
     * @throws DocumentException
     */
    static public void save_jfreechart_as_pdf(JFreeChart chart, String filename) throws IOException, DocumentException {
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

    /**
     * Default font size based on the dimension of the chart. It will overwirte the previous setting
     * @param chart
     * @return updated chart
     */
    public static JFreeChart apply_scalable_fontsize(JFreeChart chart) {
        XYPlot plot = (XYPlot) chart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();

        float baseFontSize = basicFontSize();

        Font titleFont     = new Font("Times New Roman", Font.BOLD,  (int)(baseFontSize * 1.4));
        Font labelFont     = new Font("Times New Roman", Font.BOLD,  (int)(baseFontSize * 1.2));
        Font tickFont      = new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize));
        //Font statFont      = new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize * 0.9));

        chart.getTitle().setFont(titleFont);
        xAxis.setLabelFont(labelFont);
        xAxis.setTickLabelFont(tickFont);
        yAxis.setLabelFont(labelFont);
        yAxis.setTickLabelFont(tickFont);
        chart.getLegend().setItemFont(new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize * 0.9)));

        return chart;
    }



}