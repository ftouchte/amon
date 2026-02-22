package io.github.ftouchte;
/* 
 * Code inspired by https://code.jlab.org/hallb/clas12/coatjava/jnp/-/blob/master/jnp-hipo4/src/main/java/org/jlab/jnp/hipo4/utils/HipoUtilities.java?ref_type=heads
 *
 * A routine to select the specified events from a HIPO file
 * Apply your own cuts: 
 * 		- search for : selection cuts
 * 
 * Run on an IDE or
 * Run : java -cp coat-libs-13.7.0.jar EventFilter.java
 */

import java.util.ArrayList;
import java.util.List;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriterSorted;
import org.jlab.jnp.utils.benchmark.ProgressPrintout;
import java.nio.file.Files;
import java.nio.file.Path;
import java.io.IOException;


public class EventFilter {
	public static void main(String[] args) {
		// Read filename from arguments
		// String filename = args[0];
		// String output = args[1];
		
		String filename = "/home/touchte-codjo/Desktop/hipofiles/kalman-filter/rec-data-r22712-v86.hipo";
		//String filename = "/home/touchte-codjo/Desktop/hipofiles/simulation/kalmanFilterTest/rec-simu-deuteron-v80.hipo";
		String output = "out.hipo";

		Path path = Path.of(output);

		try {
			Files.deleteIfExists(path);
		} catch (IOException e) {
			e.printStackTrace();
		}


		// Initialise HIPO reader for the file
		HipoReader reader = new HipoReader();
		reader.open(filename);

		// Read list of schemas from reader
		// If you want, you can only select the desired schemas
		// e.g create an empty List<Schema> and add
		// list.add(new Schema("AHDC::kftrack", 23000, 26))
		// refers to the bank definitions in coatjava
		SchemaFactory factory = reader.getSchemaFactory();
		List<Schema>   schemaList   = factory.getSchemaList(); // here we read all schemas

		// Create new (empty) banks from all choosed schemas
		// The will be loaded/filled during the event process
	    List<Bank>     schemaBanks  = new ArrayList<Bank>();
		for (Schema schema : schemaList){
			schemaBanks.add(new Bank(schema));
		}

		// Initialise a writer with the selected schemas
		// used to create a new HIPO file
		HipoWriterSorted writer = new HipoWriterSorted();
		for(Schema schema : schemaList){
				writer.getSchemaFactory().addSchema(schema);
		}
		// Open the output file
		writer.open(output);

		// ------------------------------------------
		// Start analysis : event selection
		// You can add your own cuts
		// ------------------------------------------
		Event   inputEvent = new Event();
		Event   outEvent   = new Event();
		ProgressPrintout progress = new ProgressPrintout();
		int counter = 0;
		int ntracks = 0;
		
		// Bank to be loaded
		Bank runConfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));
		Bank trackBank = new Bank(reader.getSchemaFactory().getSchema("AHDC::kftrack"));
		Bank matchingBank = new Bank(reader.getSchemaFactory().getSchema("ALERT::ai:projections")); // matching of AHDC and ATOF
		Bank atofHitBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::hits"));
		Bank mcBank = new Bank(reader.getSchemaFactory().getSchema("MC::Particle"));

		// Loop over events
		while(reader.hasNext()==true){ // just check if we have a next event but do not load it
			// reset outEvent
			outEvent.reset();

			// load the event
			reader.nextEvent(inputEvent); 

			// load banks to apply some cuts
			inputEvent.read(runConfigBank);
			inputEvent.read(matchingBank);
			inputEvent.read(trackBank);
			inputEvent.read(atofHitBank);
			inputEvent.read(mcBank);

			// ---------------------------
			// selection cuts
			// ---------------------------
			// // example
			// int event_number = runConfigBank.getInt("event", 0);
			// if (!listOfEvents.contains(event_number)) { 
			// 	continue; // ignore this event
			// }
			// System.out.println("event number : " + event_number );
			if (trackBank.getRows() > 0) ntracks++; // only count one track
			boolean criteria = false;
			int trackid = -1;
			int atofid = -1;
			for (int row = 0; row < matchingBank.getRows(); row++) {
				trackid = matchingBank.getInt("trackid", row);
				atofid = matchingBank.getInt("matched_atof_hit_id", row);
				if (trackid > 0 && atofid > 0) {
					criteria = true;
					break; // break the loop, we will keep this event
				}
			}
			/*if (trackid > 0 && atofid > 0) {
				int sector = -1;
				int layer = -1;
				
				// read sector/layer for this wedge
				for (int hitRow = 0; hitRow < atofHitBank.getRows(); hitRow++) {
					if (atofHitBank.getShort("id", hitRow) == atofid) {
						sector = atofHitBank.getInt("sector", hitRow);
						layer = atofHitBank.getInt("layer", hitRow);
					}
				}
				// find the corresponding bar
				if (sector > 0 && layer > 0) {
					for (int hitRow = 0; hitRow < atofHitBank.getRows(); hitRow++) {
						if (atofHitBank.getInt("sector", hitRow) == sector && atofHitBank.getInt("layer", hitRow) == layer && atofHitBank.getInt("component", hitRow) == 10) {
							criteria = true;
						}
					}
				}

				// track study
				// for (int row = 0; row < trackBank.getRows(); row++) {
				// 	if (trackid == trackBank.getInt("trackid", row)) {
				// 		double vz = 0.1*trackBank.getFloat("z", row); // cm
				// 		// get deuteron vertex
				// 		double mc_vz = mcBank.getFloat("vz", 1);
				// 		if (Math.abs(vz-mc_vz) > 1) {
				// 			criteria = criteria && false; // ignore this event
				// 		}
				// 	}
				// }
			}*/

			if (counter > 0) continue; // select only one event
			
			//if (trackBank.getInt("n_hits", 0) < 9) continue;
			if (!criteria) continue; // if the criteria is not good, ignore this event
			int event_number = runConfigBank.getInt("event", 0);
			//if (event_number != 117 || counter > 0) continue;
			System.out.println("> event : " + event_number);
			counter++;
			// ---------------------------
			// Main routine
			// ---------------------------

			// loop over all existing banks in the input file (reader)
			for(int b = 0; b < schemaBanks.size(); b++){
				// load the corresponding bank as did for runConfigBank
				inputEvent.read(schemaBanks.get(b));
				// write this bank in the output file it has at least one row
				if (schemaBanks.get(b).getRows() > 0) {
					outEvent.write(schemaBanks.get(b));
				}
			}
			// add this event to the writer
			int tag = inputEvent.getEventTag();
			outEvent.setEventTag(tag);
			writer.addEvent(outEvent,outEvent.getEventTag());
			// arbitrary
			progress.updateStatus();
		}
		System.out.println("Number of events : " + reader.getEventCount());
		System.out.println("Number of tracks : " + ntracks);
		System.out.println("Number of matches: " + counter + "==> " + 100.0*counter/ntracks + " %%");
		
		// ----------------
		// Save file
		// ----------------
		writer.close();
	}
}
