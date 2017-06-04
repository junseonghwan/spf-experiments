package pmcmc;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import briefj.BriefIO;

public class LogZProcessor<P extends ModelParameters> {

	private List<String> logZs = new ArrayList<>();
	
	public void output(File file) {
		PrintWriter writer = BriefIO.output(file);
		for (String line : logZs)
		{
			writer.println(line);
		}
		writer.close();
	}

	public void process(P p, double logZ) {
		String line = p.asCommaSeparatedLine() + ", " + logZ;
		logZs.add(line);
	}

}
