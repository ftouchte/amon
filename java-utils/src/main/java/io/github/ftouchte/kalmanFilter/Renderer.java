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
import org.jfree.chart.axis.AxisLabelLocation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.ColumnArrangement;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.chart.title.CompositeTitle;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.data.xy.XYDataset;
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

import io.github.ftouchte.kalmanFilter.Renderer.RendererOutputType;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import java.util.ArrayList;
import java.util.Arrays;
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
        double[] counts = toDoubleArray(h.getData());

        //double count_min = Arrays.stream(counts).min().getAsDouble();
        double count_min = 0;

        // Dataset
        XYSeries series = new XYSeries(h.getName());

        series.add(xLimits[0], count_min);
        for (int i = 0; i < xValues.length; i++) {
            double lower = xLimits[i];
            double upper = xLimits[i+1];
            series.add(lower, counts[i]);
            series.add(upper, counts[i]);
        }
        series.add(xLimits[xLimits.length-1], count_min);

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
     * The user can take this piece of code as example to custimize the chart
     * @param chart working chart
     * @return updated chart
     */
    public static JFreeChart apply_default_chart_rendering(JFreeChart chart) {

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
        xAxis.setTickLabelPaint(Color.BLACK);
        xAxis.setLabelPaint(Color.BLACK);
        
        yAxis.setTickUnit(new NumberTickUnit((int) (yMax - yMin) / nTicksY));
        yAxis.setAutoRangeIncludesZero(true);  // inclure 0 dans l'auto range
        yAxis.setTickLabelPaint(Color.BLACK);
        yAxis.setLabelPaint(Color.BLACK);  
        
        // control the size
        chart = apply_scalable_fontsize(chart);

        return chart;
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

        float baseFontSize = baseFontSize();

        Font titleFont     = new Font("Times New Roman", Font.BOLD,  (int)(baseFontSize * 1.4));
        Font labelFont     = new Font("Times New Roman", Font.BOLD,  (int)(baseFontSize * 1.2));
        Font tickFont      = new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize));
        //Font statFont      = new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize * 0.9));

        chart.getTitle().setFont(titleFont);
        xAxis.setLabelFont(labelFont);
        xAxis.setTickLabelFont(tickFont);
        yAxis.setLabelFont(labelFont);
        yAxis.setTickLabelFont(tickFont);
        if (chart.getLegend() != null)
            chart.getLegend().setItemFont(new Font("Times New Roman", Font.PLAIN, (int)(baseFontSize * 0.9)));

        return chart;
    }

    /**
     * To be improved !!! Do not work well yet
     * @param chart
     * @return
     */
    public static JFreeChart apply_scalable_ticksize(JFreeChart chart) {
        XYPlot plot = (XYPlot) chart.getPlot();
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();

        float baseFontSize = baseFontSize();

        float tickSize      = baseFontSize * 0.8f;
        float minorTickSize = baseFontSize * 0.4f;

        xAxis.setTickMarkInsideLength(tickSize);
        xAxis.setTickMarkOutsideLength(0.0f);
        xAxis.setTickMarkStroke(new BasicStroke(1.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER));
        xAxis.setMinorTickMarkInsideLength(minorTickSize);
        xAxis.setMinorTickMarkOutsideLength(0.0f);
        xAxis.setAxisLineStroke(new BasicStroke(1.5f));

        yAxis.setTickMarkInsideLength(tickSize);
        yAxis.setTickMarkOutsideLength(0.0f);
        yAxis.setTickMarkStroke(new BasicStroke(1.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER));
        yAxis.setMinorTickMarkInsideLength(minorTickSize);
        yAxis.setMinorTickMarkOutsideLength(0.0f);
        yAxis.setAxisLineStroke(new BasicStroke(1.5f));
        

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

        TextTitle titleLine   = new TextTitle("Stats", new Font("Monospaced", Font.BOLD, (int) (1.1*baseFontSize())));
        //TextTitle entriesLine = new TextTitle("Entries : " + h.getEntries(),               new Font("Monospaced", Font.PLAIN, (int) baseFontSize()));
        TextTitle meanLine    = new TextTitle(String.format("Mean    : %.4f", h.getMean()), new Font("Monospaced", Font.PLAIN, (int) baseFontSize()));
        TextTitle stdLine     = new TextTitle(String.format("Std     : %.4f", h.getRMS()),new Font("Monospaced", Font.PLAIN, (int) baseFontSize()));

        //System.out.println("h entries : " + h.getEntries());

        titleLine.setPosition(RectangleEdge.TOP);
        //entriesLine.setPosition(RectangleEdge.TOP);
        meanLine.setPosition(RectangleEdge.TOP);
        stdLine.setPosition(RectangleEdge.TOP);

        container.add(titleLine);
        //container.add(entriesLine);
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

    /** Color palette */
    public static Color[] colorPalette = {Color.BLACK, Color.RED, Color.GREEN, Color.BLUE, Color.MAGENTA, Color.CYAN};

    public static void restore_default_colorPalette() {
        colorPalette = new Color[] {Color.BLACK, Color.RED, Color.GREEN, Color.BLUE, Color.MAGENTA, Color.CYAN};
    }

    /**
     * Return the n-th color from the color palette. If n is bigger than colorPalette.length,
     * it returns the color associated with n modulo colorPalette.length
     * @param n-th 
     * @return Color
     */
    public static Color pickColor(int n) {
        return colorPalette[n % colorPalette.length];
    }

    /**
     * Generate n colors
     * @param n number of color
     */
    public static void generateDefaultColorPalette(int n) {
        Color[] colors = new Color[n];
        for (int i = 0; i < n; i++) {
            float hue = (float) i / n; // 0.0 (rouge) → 1.0 (rouge again)
            colors[i] = Color.getHSBColor(hue, 1.0f, 0.9f);
        }
        colorPalette = colors;
    }

    /**
     * Generate n colors
     * Jet (bleu → cyan → vert → jaune → rouge) — le classique ROOT/matplotlib
     * @param n number of color
     */
    public static void generateJetPalette(int n) {
        Color[] colors = new Color[n];
        for (int i = 0; i < n; i++) {
            float t = (float) i / (n - 1);
            float hue = (1.0f - t) * 0.66f; // 0.66 (bleu) → 0.0 (rouge)
            colors[i] = Color.getHSBColor(hue, 1.0f, 0.9f);
        }
        colorPalette = colors;
    }

    /**
     * Generate n colors
     * Turbo (meilleur que Jet pour la lisibilité) — évite le jaune trop clair
     * @param n number of color
     */
    public static void generateTurboPalette(int n) {
        Color[] colors = new Color[n];
        for (int i = 0; i < n; i++) {
            float t = (float) i / (n - 1);
            float hue        = (1.0f - t) * 0.75f; // violet → rouge
            float saturation = 0.9f;
            float brightness = 0.5f + 0.4f * (float) Math.sin(Math.PI * t); // pic au milieu
            colors[i] = Color.getHSBColor(hue, saturation, brightness);
        }
        colorPalette = colors;
    }

    /**
     * Generate n colors
     * Viridis (le meilleur pour publications scientifiques) — violet → bleu → vert → jaune
     * @param n number of color
     */
    public static void generateViridisPalette(int n) {
        // Points de contrôle RGB du vrai Viridis
        int[][] ctrl = {
            { 68,   1,  84},  // violet foncé
            { 72,  40, 120},  // violet
            { 59,  82, 139},  // bleu
            { 33, 145, 140},  // cyan
            { 50, 180,  90},  // vert moyen  ← remplacer par jaune-vert
            {180, 220,  50},  // jaune-vert  ← casse la répétition du vert
            {253, 231,  37}   // jaune
        };
        Color[] colors = new Color[n];
        for (int i = 0; i < n; i++) {
            float t    = (float) i / (n - 1) * (ctrl.length - 1);
            int   idx  = (int) t;
            float frac = t - idx;
            if (idx >= ctrl.length - 1) idx = ctrl.length - 2;
            int r = (int)(ctrl[idx][0] + frac * (ctrl[idx+1][0] - ctrl[idx][0]));
            int g = (int)(ctrl[idx][1] + frac * (ctrl[idx+1][1] - ctrl[idx][1]));
            int b = (int)(ctrl[idx][2] + frac * (ctrl[idx+1][2] - ctrl[idx][2]));
            colors[i] = new Color(r, g, b);
        }
        colorPalette = colors;
    }

    /**
     * Access basi font size scaled with the current dimensions of the charts
     * @return
     */
    public static float baseFontSize() {
        return (float) Math.sqrt(width * height) / 45.0f;
    }

    /**
     * Draw one or more histogram in a chart
     * 
     * The color lines of the histograms are based on the current colorPalette
     * 
     * See also : {@link #generateDefaultColorPalette(int)}, {@link #generateJetPalette(int)}, etc...
     * @param histos
     */
    public static JFreeChart create_overlayed_histograms(ArrayList<H1F> histos) {
        /// --- create dataset
        if (histos.isEmpty() || histos == null) return null;
        XYSeriesCollection dataset = create_dataset(histos);
        
        /// --- Create XYLineChart
        JFreeChart chart = ChartFactory.createXYLineChart(histos.get(0).getTitle(), histos.get(0).getTitleX(), histos.get(0).getTitleY(), dataset, PlotOrientation.VERTICAL, false, true, false);

        /// --- Apply default aspect to the chart
        chart = apply_default_chart_rendering(chart);

        /// --- Render plot as STEP
        XYStepRenderer renderer = new XYStepRenderer();
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setRenderer(renderer);

        /// --- Style style and color
        for (int i = 0; i < histos.size(); i++) {
            renderer.setSeriesPaint(i, pickColor(i));
            renderer.setSeriesStroke(i, null);
        }
        if (histos.size() == 1) {
            renderer.setSeriesPaint(0, Color.BLUE);
            renderer.setSeriesStroke(0, null);
            // add annotation
            plot.addAnnotation(create_stat_box(histos.get(0)));
        }
        return chart;
    }

    /**
     * Superimpose one or more histograms and save them in the specified format.
     * 
     * See also {@link {@link #create_overlayed_histograms(ArrayList)}} to customise the rendering
     * @param histos
     * @param filename
     * @param outType
     * @throws IOException
     * @throws DocumentException
     */
    public static void save_overlayed_histograms(ArrayList<H1F> histos, String filename, RendererOutputType outType) throws IOException, DocumentException {
        JFreeChart chart = create_overlayed_histograms(histos);
        Renderer.save_jfreechart(chart, filename, outType);
    }

    /**
     * Draw one histogram in a chart and save it in the specified format.
     * 
     * See also {@link {@link #create_overlayed_histograms(ArrayList)}} to customise the rendering
     * 
     * @param histo
     * @param filename
     * @param outType PNG or PDF
     * @throws IOException
     * @throws DocumentException
     */
    public static void save_histogram(H1F histo, String filename, RendererOutputType outType) throws IOException, DocumentException {
        ArrayList<H1F> list = new ArrayList<>();
        list.add(histo);
        JFreeChart chart = create_overlayed_histograms(list);
        Renderer.save_jfreechart(chart, filename, outType);
    }

    /**
     * Used by {@link #create_histogram_evolution_with_parameter(ArrayList, H1F, Color, String, String, double[])}
     * 
     * The color lines of the histograms are based on the current colorPalette
     * 
     * See also : {@link #generateDefaultColorPalette(int)}, {@link #generateJetPalette(int)}, etc...
     * 
     * @param chart initial chart
     * @param barName name of the color bar
     * @param values values to be displayed
     * @return updated chart
     */
    public static JFreeChart display_color_bar(JFreeChart chart, String barName, double[] values) {
        //XYPlot plot = (XYPlot) chart.getPlot();
        //XYSeriesCollection dataset = (XYSeriesCollection) plot.getDataset();
        int n = values.length;

        // Crée une échelle de couleurs mappée sur les indices 0 à n
        // (0 = première couleur, n = dernière) avec blanc comme couleur par défaut
        // paintScale et scaleAxis ont le même range afin de centrer les labels
        LookupPaintScale paintScale = new LookupPaintScale(-0.5, n-0.5, Color.WHITE);

        // Associe chaque index i à une couleur pickColor(i)
        // ex: i=0 → bleu, i=1 → rouge, i=2 → vert
        for (int i = 0; i < n; i++) {
            paintScale.add(i-0.5, pickColor(i)); // couleur i commence à i-0.5
        }

        // Axe numérique pour la color bar
        NumberAxis scaleAxis = new NumberAxis(barName);

        // Centrer les labels par rapport aux couleurs
        // Décaler les ticks de 0.5 pour qu'ils soient au centre de chaque couleur
        scaleAxis.setRange(-0.5, n - 0.5); // décale le range de 0.5
        // Un tick tous les 1 unité → un tick par couleur
        scaleAxis.setTickUnit(new NumberTickUnit(1));
        scaleAxis.setTickLabelFont(new Font("Times New Roman", Font.PLAIN, (int) (baseFontSize() * 0.8f)));
        scaleAxis.setLabelFont(new Font("Times New Roman", Font.BOLD,  (int) (baseFontSize() * 0.9f)));

        // Remplace le format numérique par défaut par un format custom
        // qui convertit l'index → valeur réelle à afficher
        scaleAxis.setNumberFormatOverride(new java.text.NumberFormat() {

            // Appelé pour chaque tick : reçoit l'index (0.0, 1.0, 2.0...)
            // et retourne la vraie valeur à afficher (0.5, 1.0, 2.0...)
            @Override
            public StringBuffer format(double number, StringBuffer toAppendTo, java.text.FieldPosition pos) {
                int idx = (int) Math.round(number); // 0.0 → 0, 1.0 → 1, etc.
                if (idx >= 0 && idx < values.length) {
                    toAppendTo.append(values[idx]); // affiche values[0]=0.5, values[1]=1.0, etc.
                }
                return toAppendTo;
            }

            // Même chose pour les entiers (non utilisé ici mais obligatoire)
            @Override
            public StringBuffer format(long number, StringBuffer toAppendTo, java.text.FieldPosition pos) {
                return toAppendTo;
            }

            // Parse string → number (non utilisé ici mais obligatoire)
            @Override
            public Number parse(String source, java.text.ParsePosition parsePosition) {
                return null;
            }
        });

        PaintScaleLegend colorBar = new PaintScaleLegend(paintScale, scaleAxis);
        colorBar.setPosition(RectangleEdge.LEFT);  // à gauche
        colorBar.setMargin(new RectangleInsets(10, 5, 10, 5));
        //colorBar.setMargin(new RectangleInsets(0, 5, 0, 5)); 
        colorBar.setStripWidth(20);               // largeur de la barre
        colorBar.setAxisLocation(AxisLocation.BOTTOM_OR_LEFT);

        chart.addSubtitle(colorBar); // ajouter au chart

        return chart;

    }

    /**
     * Draw one or more histograms in a chart. Assuming the histogram vary with a parameter, this routine draw a color bar
     * @param histos list of histograms
     * @param h_ref reference histogram, can be null
     * @param color_ref color of the reference histogram, by default it is BLACK
     * @param title title of the figure
     * @param parameterName name of the parameter
     * @param values the parameter values associated to each histograms
     */
    public static JFreeChart create_histogram_evolution_with_parameter(ArrayList<H1F> histos, H1F h_ref, Color color_ref, String title, String parameterName, double[] values) {
        /// --- create dataset
        if (histos.isEmpty() || histos == null) return null;
        ArrayList<H1F> histos_and_ref = new ArrayList<>();
        for (H1F h : histos) {
            histos_and_ref.add(h);
        }
        if (h_ref != null){
            histos_and_ref.add(h_ref);
        }
        XYSeriesCollection dataset = create_dataset(histos_and_ref);

        /// --- Create XYLineChart
        JFreeChart chart = ChartFactory.createXYLineChart(title, histos.get(0).getTitleX(), histos.get(0).getTitleY(), dataset, PlotOrientation.VERTICAL, false, true, false);

        /// --- Apply default aspect to the chart
        chart = apply_default_chart_rendering(chart);

        /// --- Render plot as STEP
        XYStepRenderer renderer = new XYStepRenderer();
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setRenderer(renderer);

        /// --- Style style and color
        for (int i = 0; i < histos.size(); i++) {
            renderer.setSeriesPaint(i, pickColor(i)); // line color
            renderer.setSeriesStroke(i, null); // lien style
        }
        if (h_ref != null) {
            renderer.setSeriesPaint(histos.size(), color_ref != null ? color_ref : Color.BLACK);
            renderer.setSeriesStroke(histos.size(), null);
        }
        if (histos.size() == 1) {
            renderer.setSeriesPaint(0, Color.BLUE);
            renderer.setSeriesStroke(0, null);
            // add annotation
            plot.addAnnotation(create_stat_box(histos.get(0)));
            if (chart.getLegend() != null)
                chart.getLegend().setVisible(false);
        } else {
            chart = display_color_bar(chart, parameterName, values);
        }
        return chart;
    }

    /**
     * Save jfreechart as PDF or PNG
     * @param chart
     * @param filename
     * @param outType
     * @throws IOException
     * @throws DocumentException
     */
    public static void save_jfreechart(JFreeChart chart, String filename, RendererOutputType outType) throws IOException, DocumentException {
        if (outType == RendererOutputType.PNG) {
            save_jfreechart_as_png(chart, sanitizeFilenameForPNG(filename));
        }
        else if (outType == RendererOutputType.PDF) {
            save_jfreechart_as_pdf(chart, sanitizeFilenameForPDF(filename));
        }
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

    public static JFreeChart create_chart_for_graph(XYDataset dataset, String title, String xtitle, String ytitle) {
        JFreeChart chart = ChartFactory.createXYLineChart(title, xtitle, ytitle, dataset, PlotOrientation.VERTICAL, true, true, false);
        XYPlot plot = (XYPlot)chart.getPlot();
        chart = apply_default_chart_rendering(chart);
        chart = apply_scalable_fontsize(chart);

        /// Defalt rendering
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        int num_series = plot.getSeriesCount();
        for (int i = 0; i < num_series; i++) {
            renderer.setSeriesPaint(i, pickColor(i));
        }
        // renderer.setSeriesLinesVisible(0, false);
        // renderer.setSeriesShapesVisible(1, false);
        plot.setRenderer(renderer);
      return chart;
    }

    



}