package io.github.ftouchte;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.base.Component;
import org.jlab.geom.detector.alert.ATOF.AlertTOFDetector;
import org.jlab.geom.detector.alert.ATOF.AlertTOFFactory;
import org.jlab.geom.prim.Point3D;
/**
 * Code inspired by atof_factory.groovy in clas12Tags
 * 
 */
public class BuildATOF {
    public static void main(String[] args) {
        BuildATOF obj = new BuildATOF();
        obj.run();
    }

    AlertTOFDetector detector = null;

    public BuildATOF () {
        detector = (new AlertTOFFactory()).createDetectorCLAS(new DatabaseConstantProvider());
    }

    public void run () {
        Path path = Path.of("ATOF_geometry.txt");
        BufferedWriter writer = null;
        try {
            writer = Files.newBufferedWriter(path, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
            // // header
            writer.write("sector,superlayer,layer,component,z_top,z_bottom,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7");
            writer.newLine();

            for (int s = 0; s < detector.getNumSectors(); s++) {
                for (int sl = 0; sl < detector.getSector(s).getNumSuperlayers(); sl++) {
                    for (int l = 0; l < detector.getSector(s).getSuperlayer(sl).getNumLayers(); l++) {
                        if (sl == 0) {
                            Component comp = detector.getSector(s).getSuperlayer(sl).getLayer(l).getComponent(10);
                            writer.write(get_string(comp, s, sl, l, 10));
                            writer.newLine();
                        } else {
                            for (int c = 0; c < detector.getSector(s).getSuperlayer(sl).getLayer(l).getNumComponents(); c++) {
                                Component comp = detector.getSector(s).getSuperlayer(sl).getLayer(l).getComponent(c);
                                writer.write(get_string(comp, s, sl, l,c));
                                writer.newLine();
                            }
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

    String get_string(Component comp, int s, int sl, int l, int c) {
        // reading top face vertices of ATOF cell and storing their x,y coordinates
        Point3D p0 = comp.getVolumePoint(0);
        // double top_x_0 = p0.x();
        // double top_y_0 = p0.y();
        Point3D p1 = comp.getVolumePoint(1);
        // double top_x_1 = p1.x();
        // double top_y_1 = p1.y();
        Point3D p2 = comp.getVolumePoint(2);
        // double top_x_2 = p2.x();
        // double top_y_2 = p2.y();
        Point3D p3 = comp.getVolumePoint(3);
        // double top_x_3 = p3.x();
        // double top_y_3 = p3.y();
        // reading bottom face vertices of ATOF cell and storing their x,y coordinates
        Point3D p4 = comp.getVolumePoint(4);
        // double bottom_x_4 = p4.x();
        // double bottom_y_4 = p4.y();
        Point3D p5 = comp.getVolumePoint(5);
        // double bottom_x_5 = p5.x();
        // double bottom_y_5 = p5.y();
        Point3D p6 = comp.getVolumePoint(6);
        // double bottom_x_6 = p6.x();
        // double bottom_y_6 = p6.y();
        Point3D p7 = comp.getVolumePoint(7);
        // double bottom_x_7 = p7.x();
        // double bottom_y_7 = p7.y();
        // z_top face, z_bottom face, x, y of points 0 to 7
        return String.format("%2d  %2d  %2d  %2d % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f % 9.4f", s, sl, l, c, p0.z(), p4.z(), p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), p4.x(), p4.y(), p5.x(), p5.y(), p6.x(), p6.y(), p7.x(), p7.y());
    }
}
