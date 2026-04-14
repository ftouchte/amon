package io.github.ftouchte.filtering;

import org.jlab.io.base.DataEvent;

import io.github.ftouchte.utils.ParticleRow;
import io.github.ftouchte.utils.Units;

import org.jlab.io.base.DataBank;
import java.util.ArrayList;

/**
 * Inspired by amon/analysis/kalman-filter/kfmon_data.cpp
 * 
 * For now: the study is made on a deuteron target
 */
public class AlertElasticAnalyser {
    
    // Constants
    double beam_energy = 2.23951 * Units.GeV; // incident energy if the electron, GeV
    double electron_mass = 0.511e-3 * Units.GeV; // energy mass of electron, GeV
    double proton_mass = 938.272e-3 * Units.GeV; // energy mass of proton, GeV
    double helium_mass = 3.73 * Units.GeV; // energy mass of Helium-4, GeV
    double deuteron_mass = 1.875 * Units.GeV; // energy mass of Deuterium, GeV

    // Elatsics cuts
    double delta_phi = 20;
    double W2_min = 3.5 * Units.GeV; // GeV
    double W2_max = 3.8 * Units.GeV; // GeV
    int nhits_min = 7;

    // Particles
    ParticleRow electron;
    ParticleRow ahdc_track;

    /**
     * Check if this event containts an elastic couple : electron + track
     * 
     * If true, one can retrive the rows of the electron and of the track using {@link #getElectron()} and {@link #getAhdcTrack()}
     * 
     * @param event
     */
    public boolean IsElastic(DataEvent event) {
        if (!event.hasBank("REC::Particle") && !event.hasBank("AHDC::kftrack")) {
            return false;
        }
        
        DataBank recBank = event.getBank("REC::Particle");
        DataBank trackBank = event.getBank("AHDC::kftrack");

        electron = null;
        ahdc_track = null;

        
        for (int row = 0; row < recBank.rows(); row++) {
            // Select trigger electrons
            if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0) {
                // compute kinematic variables
                double px = recBank.getFloat("px", row);
                double py = recBank.getFloat("py", row);
                double pz = recBank.getFloat("pz", row);
                ParticleRow electron_candidate = new ParticleRow(px * Units.GeV, py * Units.GeV, pz * Units.GeV);
                electron_candidate.SetBankRow(row);
                // physics kinematics
                double scattered_beam_energy = Math.sqrt(Math.pow(electron_candidate.p(Units.GeV),2) + Math.pow(electron_mass,2));
                double nu = beam_energy - scattered_beam_energy;
                double Q2 = 4*beam_energy*scattered_beam_energy*Math.pow(Math.sin(electron_candidate.theta(Units.rad)/2),2);
                double W2 = Math.pow(deuteron_mass,2) + 2*deuteron_mass*nu - Q2;
                // W2 cut to retrieve the target mass
                if (W2 > W2_min && W2 < W2_max) {
                    electron = electron_candidate;
                    // Select the best AHDC track
                    ArrayList<ParticleRow> all_tracks = new ArrayList<>();
                    for (int trow = 0; trow < trackBank.rows(); trow++) {
                        double px1 = trackBank.getFloat("px", trow);
                        double py1 = trackBank.getFloat("py", trow);
                        double pz1 = trackBank.getFloat("pz", trow);
                        ParticleRow track_candidate = new ParticleRow(px1 * Units.MeV, py1 * Units.MeV, pz1 * Units.MeV);
                        track_candidate.SetBankRow(trow);
                        if (trackBank.getInt("n_hits", trow) >= nhits_min) {
                            all_tracks.add(track_candidate);
                        }
                    }
                    // Sort track accordingly to the smallest delta phi
                    all_tracks.sort((t1, t2) -> {
                        double dphi1 = deltaPhi(electron.phi(Units.deg), t1.phi(Units.deg));
                        double dphi2 = deltaPhi(electron.phi(Units.deg), t2.phi(Units.deg));
                        return Double.compare(dphi1, dphi2);
                    });
                    // Select the best track
                    if (all_tracks.size() > 0) {
                        ParticleRow best_track = all_tracks.get(0);
                        if (deltaPhi(electron.phi(Units.deg), best_track.phi(Units.deg)) < delta_phi) {
                            ahdc_track = best_track;
                            return true;
                        }
                    }
                    
                }
            }
        }
        
        if (electron != null && ahdc_track != null) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Check if this event containts an elastic electron. 
     * 
     * If true, one can retrive the rows of the electron {@link #getElectron()}. The theoretical kinematics of the track can be retrieved using {@link #getAhdcTrack()}
     * 
     * @param event
     */
    public boolean hasElasticElectron(DataEvent event) {
        if (!event.hasBank("REC::Particle") && !event.hasBank("AHDC::kftrack")) {
            return false;
        }
        
        DataBank recBank = event.getBank("REC::Particle");
        //DataBank trackBank = event.getBank("AHDC::kftrack");

        electron = null;
        ahdc_track = null;

        
        for (int row = 0; row < recBank.rows(); row++) {
            // Select trigger electrons
            if (recBank.getInt("pid", row) == 11 && recBank.getShort("status", row) < 0) {
                // compute kinematic variables
                double px = recBank.getFloat("px", row);
                double py = recBank.getFloat("py", row);
                double pz = recBank.getFloat("pz", row);
                ParticleRow electron_candidate = new ParticleRow(px * Units.GeV, py * Units.GeV, pz * Units.GeV);
                // physics kinematics
                double scattered_beam_energy = Math.sqrt(Math.pow(electron_candidate.p(Units.GeV),2) + Math.pow(electron_mass,2));
                double nu = beam_energy - scattered_beam_energy;
                double Q2 = 4*beam_energy*scattered_beam_energy*Math.pow(Math.sin(electron_candidate.theta(Units.rad)/2),2);
                double W2 = Math.pow(deuteron_mass,2) + 2*deuteron_mass*nu - Q2;
                // W2 cut to retrieve the target mass
                if (W2 > W2_min && W2 < W2_max) {
                    electron = electron_candidate;
                    // compute theoretical ahdc track
                    double px1 = -scattered_beam_energy*Math.sin(electron.theta(Units.rad))*Math.cos(electron.phi(Units.rad));
                    double py1 = -scattered_beam_energy*Math.sin(electron.theta(Units.rad))*Math.sin(electron.phi(Units.rad));
                    double pz1 = beam_energy - scattered_beam_energy*Math.cos(electron.theta(Units.rad));
                    ahdc_track = new ParticleRow(px1, py1, pz1);
                    return true;           
                }
            }
        }
        
        if (electron != null && ahdc_track != null) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * 
     * @param electron_phi in deg
     * @param track_phi in deg
     */
    private double deltaPhi(double electron_phi, double track_phi) {
        return Math.abs(Math.abs(electron_phi - track_phi) - 180);
    }

    /**
     * Return elactic electron. Better use after {@link #IsElastic(DataEvent)}
     */
    public ParticleRow getElectron() {
        return this.electron;
    }

    /**
     * Return elactic ahdc track. Better use after {@link #IsElastic(DataEvent)}
     */
    public ParticleRow getAhdcTrack() {
        return this.ahdc_track;
    }
}
