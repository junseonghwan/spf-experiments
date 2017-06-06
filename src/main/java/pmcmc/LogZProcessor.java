package pmcmc;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import briefj.BriefIO;

public class LogZProcessor<P extends ModelParameters> {

	private List<String> logZs = new ArrayList<>();
	private String outputPrefix;
	
	public LogZProcessor(String outputPrefix) { this.outputPrefix = outputPrefix; }

	public void output(File file) {
		PrintWriter writer = BriefIO.output(file);
		for (String line : logZs)
		{
			writer.println(line);
		}
		writer.close();
	}

	public void process(int currentIteration, P p, double logZ) {
		String line = currentIteration + ", " + p.asCommaSeparatedLine() + ", " + logZ;
		logZs.add(line);
	}
	
	public String getOutputPrefix() { return outputPrefix; }

}
