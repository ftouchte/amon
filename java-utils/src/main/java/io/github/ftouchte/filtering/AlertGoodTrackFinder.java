package io.github.ftouchte.filtering;

import java.util.ArrayList;
import org.jlab.io.base.DataEvent;

import org.jlab.io.base.DataBank;

public class AlertGoodTrackFinder implements AlertTrackSelector {
    

    public boolean hasGoodTrack(DataEvent event) {
        return method2(event);
    }


    public ArrayList<int[]> getAhdcKFTrackRows() {
        return trackRows;
    }

    private ArrayList<int[]> trackRows = new ArrayList<>();

    /**
     * Here we look at dEdx versus p and we apply some cuts
     *    - cut on p
     *    - cut on dEdx
     *    - cut on n_hits
     * @param event
     * @return
     */
    // private boolean method1(DataEvent event) {
    //     if (!event.hasBank("AHDC::kftrack")) return false;

    //     DataBank trackBank = event.getBank("AHDC::kftrack");
        
    //     trackRows.clear();

    //     for (int row = 0; row < trackBank.rows(); row++) {
    //         double px = trackBank.getFloat("px", row); // MeV
    //         double py = trackBank.getFloat("py", row); // MeV
    //         double pz = trackBank.getFloat("pz", row); // MeV
    //         double p = Math.sqrt(Math.pow(px, 2)+Math.pow(py, 2)+Math.pow(pz, 2)); // MeV
    //         double dEdx = trackBank.getFloat("dEdx", row);
    //         int nhits = trackBank.getInt("n_hits", row);
    //         double vz = trackBank.getFloat("z", row)*0.1; // cm
    //         double chi2 = trackBank.getFloat("chi2", row);

            
    //         if ((nhits >= 7 && vz > -15 && vz < 15 && chi2 > 0.5 && chi2 < 1.5)) {
    //             trackRows.add(row);
    //         }
    //     }

    //     return trackRows.size() > 0;
    // }

    /**
     * An electron with a good vz (so no FT) + a track with chi2 < 8
     * @param event
     * @return
     */
    private boolean method2(DataEvent event) {

        if (!event.hasBank("REC::Particle") && !event.hasBank("AHDC::kftrack")) {
            return false;
        }
        
        DataBank recBank = event.getBank("REC::Particle");
        DataBank trackBank = event.getBank("AHDC::kftrack");

        trackRows.clear();

        // Select electron
        for (int row = 0; row < recBank.rows(); row++) {
            if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0) {

                int status = recBank.getShort("status", row);

                // Ignore FT electrons
                if (Math.abs(status)/1000 == 1) { 
                    continue;
                }

                // Select track
                for (int i = 0; i < trackBank.rows(); i++) {
                    int nhits = trackBank.getInt("n_hits", i);
                    double chi2 = trackBank.getFloat("chi2", i);
                    
                    // Select good track
                    if (nhits < 6 || chi2 > 8) continue;
                    int[] vec = {row, i};
                    trackRows.add(vec);
                }

            }
        }

        return false;

    }

}
