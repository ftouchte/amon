package io.github.ftouchte.filtering;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriterSorted;

/**
 * Filter elastic events from file.
 * 
 * Please refer to {@link AlertElasticAnalyser#AlertElasticAnalyser}
 */
public class ElasticEventFilter {
    public static void main(String[] args) {

        if (args.length < 2) {
            System.out.println("* error : missing arguments");
            System.out.println("* require an input file and an output file");
        }

        // Read input and ouputs file from arguments
        // String inFile = "/volatile/clas12/touchte/kalman-filter/v120/rec_clas_022712.evio.00000.hipo";
		// String outFile = "/lustre24/expphy/volatile/clas12/touchte/kalman-filter/test/elastic_" + Path.of(inFile).getFileName().toString();
        String inFile = args[0];
		String outFile = args[1];

        
		Path path = Path.of(outFile);
        
        // Check that the input file exists
        if (Files.exists(Path.of(inFile))) {
            System.out.println("> Input file : " + inFile);
		} else {
            System.out.println("* " + path + " : error file not found");
            return;
        }
        
        // Check that the provided output file does not exist
        if (Files.exists(path)) {
			System.out.println(path + " : error existing file");
            return;
		}

        // Initialise HIPO reader for the file
		HipoReader reader = new HipoReader();
		reader.open(inFile);

        // Read list of schemas from reader
		// If you want, you can only select the desired schemas
		// e.g create an empty List<Schema> and add
		// list.add(new Schema("AHDC::kftrack", 23000, 26))
		// in that case, refer to the bank definitions in coatjava
		SchemaFactory factory = reader.getSchemaFactory();
		List<Schema>   schemaList   = factory.getSchemaList();

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
		writer.open(outFile);

        // ------------------------------------------
		// Start analysis : event selection
		// You can add your own cuts
		// ------------------------------------------
		Event   inEvent = new Event();
		Event   outEvent   = new Event();

        // Initialiase elastic analyser
        AlertElasticAnalyser analyser = new AlertElasticAnalyser();

        // Loop over events
        int counter = 0;
        int nevents = 0;
        while(reader.hasNext()==true){ // just check if we have a next event but do not load it
            // print
            nevents++;
            if (nevents % 10000 == 0) {
                System.out.println("  * processed events : " + nevents);
            }
            
            // reset outEvent
			outEvent.reset();

			// load the event
			reader.nextEvent(inEvent);

            // Check that the event is an elastic
            if (analyser.IsElastic(new HipoDataEvent(inEvent, factory))) {
                counter++;
                // loop over all existing banks in the input file (reader)
                for(int b = 0; b < schemaBanks.size(); b++){
                    // load the corresponding bank as did for runConfigBank
                    inEvent.read(schemaBanks.get(b));
                    // write this bank in the output file it has at least one row
                    if (schemaBanks.get(b).getRows() > 0) {
                        outEvent.write(schemaBanks.get(b));
                    }
                }
                // add this event to the writer
                int tag = inEvent.getEventTag();
                outEvent.setEventTag(tag);
                writer.addEvent(outEvent,outEvent.getEventTag());
            }

        } // end loop over events

        // ----------------
		// Save file
		// ----------------
		writer.close();
        System.out.println("> Output file : " + outFile);
        System.out.println("> Nb elastics : " + counter);
    }
}
