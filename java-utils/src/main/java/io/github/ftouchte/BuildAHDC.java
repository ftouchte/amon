package io.github.ftouchte;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCWire;
import org.jlab.geom.prim.Point3D;
/**
 * Code inspired by atof_factory.groovy in clas12Tags
 * 
 */
public class BuildAHDC {
    public static void main(String[] args) {
        BuildAHDC obj = new BuildAHDC();
        obj.run();
    }

    AlertDCDetector detector = null;

    public BuildAHDC () {
        detector = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());
    }

    public void run () {
        Path path = Path.of("AHDC_geometry.txt");
        BufferedWriter writer = null;
        try {
            writer = Files.newBufferedWriter(path, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
            // // header
            writer.write("sector,superlayer,layer,component,z_origin,z_end,x1,y1,x2,y2");
            writer.newLine();

            for (int s = 1; s <= detector.getNumSectors(); s++) {
                for (int sl = 1; sl <= detector.getSector(s).getNumSuperlayers(); sl++) {
                    for (int l = 1; l <= detector.getSector(s).getSuperlayer(sl).getNumLayers(); l++) {
                        for (int c = 1; c <= detector.getSector(s).getSuperlayer(sl).getLayer(l).getNumComponents(); c++) {
                            AlertDCWire comp = detector.getSector(s).getSuperlayer(sl).getLayer(l).getComponent(c);
                            writer.write(get_string(comp, s, sl, l,c));
                            writer.newLine();
                        } 
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("File created : " + path.toString());
        
        
    }

    String get_string(AlertDCWire comp, int s, int sl, int l, int c) {
        // reading top face vertices of ATOF cell and storing their x,y coordinates
        Point3D origin = comp.getLine().origin();
        Point3D end = comp.getLine().end();
        String str = String.format("%2d  %2d  %2d  %2d % 8.4f % 8.4f % 8.4f % 8.4f % 8.4f % 8.4f", s, sl, l, c, origin.z(), end.z(), origin.x(), origin.y(), end.x(), end.y());
        System.out.println(str);
        return str;
    }
}
