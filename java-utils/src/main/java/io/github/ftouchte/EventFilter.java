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


public class EventFilter {
	public static void main(String[] args) {
		// Read filename from arguments
		// String filename = args[0];
		// String output = args[1];
		
		String filename = "/home/touchte-codjo/Desktop/hipofiles/kalman-filter/rec-data-r22712-v86.hipo";
		String output = "out.hipo";


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
		int track_counter = 0;
		
		// Bank to be loaded
		Bank runConfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));
		Bank trackBank = new Bank(reader.getSchemaFactory().getSchema("AHDC::kftrack"));
		Bank matchingBank = new Bank(reader.getSchemaFactory().getSchema("ALERT::ai:projections")); // matching of AHDC and ATOF

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

			// ---------------------------
			// selection cuts
			// ---------------------------
			// // example
			// int event_number = runConfigBank.getInt("event", 0);
			// if (!listOfEvents.contains(event_number)) { 
			// 	continue; // ignore this event
			// }
			// System.out.println("event number : " + event_number );
			// if (counter > 0) break;
			//if (trackBank.getRows() > 0) track_counter++;
			if (trackBank.getInt("n_hits", 0) < 9) continue;
			boolean criteria = false;
			for (int row = 0; row < matchingBank.getRows(); row++) {
				int trackid = matchingBank.getInt("trackid", row);
				int atof_wedge_id = matchingBank.getInt("matched_atof_hit_id", row);
				if (trackid > 0 && atof_wedge_id > 0) {
					criteria = true;
					break; // the loop, we will keep this event
				}
			}
			if (!criteria) continue; // if the criteria is not good, ignore this event
			
			counter++;
			int event_number = runConfigBank.getInt("event", 0);
			System.out.println("> event : " + event_number);
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
		System.out.println("Number of tracks : " + track_counter);
		System.out.println("Number of matches: " + counter);
		System.out.printf ("Percentage       : %.2f %%\n", 1.0*counter/track_counter);
		
		// ----------------
		// Save file
		// ----------------
		writer.close();
	}
}
