package io.github.ftouchte.filtering;

import java.util.ArrayList;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;

public class AlertGoodTrackFinder implements AlertTrackSelector {
    

    public boolean hasGoodTrack(DataEvent event) {
        return method1(event);
    }


    public ArrayList<Integer> getAhdcKFTrackRows() {
        return trackRows;
    }

    private ArrayList<Integer> trackRows = new ArrayList<>();

    /**
     * Here we look at dEdx versus p and we apply some cuts
     *    - cut on p
     *    - cut on dEdx
     *    - cut on n_hits
     * @param event
     * @return
     */
    private boolean method1(DataEvent event) {
        if (!event.hasBank("AHDC::kftrack")) return false;

        DataBank trackBank = event.getBank("AHDC::kftrack");
        
        trackRows.clear();

        for (int row = 0; row < trackBank.rows(); row++) {
            double px = trackBank.getFloat("px", row); // MeV
            double py = trackBank.getFloat("py", row); // MeV
            double pz = trackBank.getFloat("pz", row); // MeV
            double p = Math.sqrt(Math.pow(px, 2)+Math.pow(py, 2)+Math.pow(pz, 2)); // MeV
            double dEdx = trackBank.getFloat("dEdx", row);
            int nhits = trackBank.getInt("n_hits", row);
            double vz = trackBank.getFloat("z", row)*0.1; // cm
            double chi2 = trackBank.getFloat("chi2", row);

            
            if ((nhits >= 7 && vz > -15 && vz < 15 && chi2 > 0.5 && chi2 < 1.5)) {
                trackRows.add(row);
            }
        }

        return trackRows.size() > 0;
    }
}
