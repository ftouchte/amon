package io.github.ftouchte.utils;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;

import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * Extrait un GraphErrors stocké dans un fichier .hipo (via GROOT/TDirectory)
 * et l'exporte en CSV (x, y, ex, ey) pour relecture facile en C++/ROOT.
 *
 * Compilation (avec le jar groot dans le classpath) :
 *   javac -cp groot.jar ExtractGraph.java
 * Execution :
 *   java -cp .:groot.jar ExtractGraph
 */
public class ExtractGraph {

    public static void main(String[] args) throws Exception {

        // -------- A REMPLIR --------
        String hipoFile   = "/volatile/clas12/rg-l/production/tline/pass0_v10.3_alert/timeline_web/forward/forward_electron_VZ.hipo";   // chemin du fichier .hipo
        String folderPath = "timelines/";         // chemin du dossier dans le hipo (vu dans hipo-browser)
        String graphName  = "sec1_upstream";     // nom exact du dataset (vu dans hipo-browser)
        String outputCsv  = "graph_export.csv";
        // ---------------------------

        TDirectory dir = new TDirectory();
        dir.readFile(hipoFile);

        GraphErrors graph = (GraphErrors) dir.getObject(folderPath, graphName);

        if (graph == null) {
            System.err.println("Graphe introuvable : verifie folderPath et graphName avec hipo-browser.");
            return;
        }

        int n = graph.getDataSize(0);
        System.out.println("Nombre de points : " + n);

        try (PrintWriter pw = new PrintWriter(new FileWriter(outputCsv))) {
            pw.println("x,y,ex,ey");
            for (int i = 0; i < n; i++) {
                double x  = graph.getDataX(i);
                double y  = graph.getDataY(i);
                double ex = graph.getDataEX(i);
                double ey = graph.getDataEY(i);
                pw.printf("%g,%g,%g,%g%n", x, y, ex, ey);
            }
        }

        System.out.println("Export termine -> " + outputCsv);
    }
}