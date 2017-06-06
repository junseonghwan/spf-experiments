package experiments.kitagawa;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import briefj.BriefIO;
import pmcmc.PMCMCProcessor;

public class KitagawaProcessor implements PMCMCProcessor<KitagawaParams> {

	private List<String> lines = new ArrayList<>();
	private String outputPrefix;
	public KitagawaProcessor(String outputPrefix) {
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
	public void process(KitagawaParams p) {
		lines.add(p.asCommaSeparatedLine());
	}

	@Override
	public String outputPrefix() {
		return outputPrefix;
	}

}
