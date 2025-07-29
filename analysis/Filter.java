/* 
 * Inspired by https://code.jlab.org/hallb/clas12/coatjava/jnp/-/blob/master/jnp-hipo4/src/main/java/org/jlab/jnp/hipo4/utils/HipoUtilities.java?ref_type=heads
 *
 * For the moment, the output file is out.hipo
 *
 * Run : java -cp coat-libs-13.1.0.jar Filter.java
 */

import java.io.Console;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoChain;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.hipo4.io.HipoWriterSorted;
import org.jlab.jnp.utils.benchmark.Benchmark;
import org.jlab.jnp.utils.benchmark.ProgressPrintout;
import org.jlab.jnp.utils.data.ArrayUtils;
import org.jlab.jnp.utils.options.OptionStore;
import org.jlab.jnp.utils.file.FileUtils;


public class Filter {
	public static void main(String[] args) {
		List<Integer> listOfEvents = new ArrayList<>(Arrays.asList(173000, 172985, 172981, 173023, 173004, 173038, 173049, 172927, 173057, 173091, 173101, 173041, 173058, 173059, 173042, 173105, 173126, 173129, 173131, 173123, 173111, 172922, 173116, 173150, 173134, 173110, 173195, 173190, 173185, 173167, 173201, 173194, 173248, 173252, 173094, 173271, 173240, 173277, 173280, 173292, 173324, 173293, 173331, 173339, 173294, 173377, 173389, 173361, 173349, 173445, 173262, 173431, 173442, 173447, 173472, 173525, 173553, 173550, 173511, 173506, 173523, 173492, 173586, 173596, 173626, 173627, 173656, 173673, 173612, 173693, 173713, 173719, 173593, 173684, 173752, 173784, 173786, 173769, 173792, 173751, 173761, 173823, 173762, 173819, 173827, 173822, 173771, 173851, 173733, 173857, 173833, 173875, 173843, 173890, 173716, 173903));
		System.out.println("Test: coat-libs-13.1.0.jar");
		String filename = "/home/touchte-codjo/Desktop/hipofiles/occupancy/rec_clas_022435.evio.00003.hipo";
		HipoReader reader = new HipoReader();
		reader.open(filename);
		Event event = new Event();
		// pre start (to be understood)
		SchemaFactory factory = reader.getSchemaFactory();
		HipoWriterSorted writer = new HipoWriterSorted();
		List<Schema>   schemaList   = factory.getSchemaList();
	        List<Bank>     schemaBanks  = new ArrayList<Bank>();
		for (Schema schema : schemaList){
			schemaBanks.add(new Bank(schema));
		}
	        for(Schema schema : schemaList){
                	writer.getSchemaFactory().addSchema(schema);
        	}
		// start (to be understood)
		writer.open("out.hipo");
	        reader.close();
		Event   inputEvent = new Event();
		Event   outEvent   = new Event();
		ProgressPrintout progress = new ProgressPrintout();
		// only one file
		HipoReader ir = new HipoReader();
		ir.open(filename);
		System.out.println("****>>>> openning file : " + filename);
		int counter = 0;
		System.out.println("--> number of events = " + ir.getEventCount());
		Bank runConfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));
		while(ir.hasNext()==true){
			outEvent.reset();
			ir.nextEvent(inputEvent);
			// my filter
			inputEvent.read(runConfigBank);
			int event_number = runConfigBank.getInt("event", 0);
			if (listOfEvents.contains(event_number)) { 
				System.out.println("event number : " + event_number );
				// end my filter
				int tag = inputEvent.getEventTag();
				outEvent.setEventTag(tag);
				for(int b = 0; b < schemaBanks.size(); b++){
					inputEvent.read(schemaBanks.get(b));
					if (schemaBanks.get(b).getRows()>0) {
						outEvent.write(schemaBanks.get(b));
					}
				}
				counter++;
				writer.addEvent(outEvent,outEvent.getEventTag());
				progress.updateStatus();
			}
			/*if (counter > 10) {
				writer.close();
				return;
			}*/
		}
		System.out.println("****>>>> processed events : " + counter);
		writer.close();
	}
}
