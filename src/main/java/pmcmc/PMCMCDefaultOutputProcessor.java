package pmcmc;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import briefj.BriefIO;

public class PMCMCDefaultOutputProcessor<P extends ModelParameters> implements PMCMCProcessor<P> {

	private List<String> lines = new ArrayList<>();
	private String outputPrefix;
	public PMCMCDefaultOutputProcessor(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	@Override
	public void output(File file) {
		PrintWriter writer = BriefIO.output(file);
		for (String line : lines)
		{
			writer.println(line);
		}
		writer.close();
	}

	@Override
	public void process(P p) {
		lines.add(p.asCommaSeparatedLine());
	}

	@Override
	public String outputPrefix() {
		return outputPrefix;
	}

}
