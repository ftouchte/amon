package io.github.ftouchte.filtering;

import java.util.ArrayList;

import org.jlab.io.base.DataEvent;

public interface AlertTrackSelector {
    
    /**
     * 
     * @param event Data event to be processed
     * @return true, if this event contains good tracks
     */
    public boolean hasGoodTrack(DataEvent event);

    /**
     * Should ideally be used just after {@link #hasGoodTrack(DataEvent)}. If not, should return an empty list.
     * @return the couple of rows for (REC::Particle ; AHDC::kftrack)  identified as good events. If electron is associated return a negative number, e.g -1
     */
    public ArrayList<int[]> getAhdcKFTrackRows();

}
