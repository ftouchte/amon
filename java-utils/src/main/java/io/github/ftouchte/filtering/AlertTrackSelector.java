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
     * Should ideally be used jsut after {@link #hasGoodTrack(DataEvent)}. If not, should return an empty list.
     * @return the rows, in AHDC::kftrack, of all tracks considered as good ones.
     */
    public ArrayList<Integer> getAhdcKFTrackRows();

}
