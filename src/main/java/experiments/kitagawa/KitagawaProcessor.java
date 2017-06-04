package experiments.kitagawa;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import briefj.BriefIO;
import pmcmc.PMCMCProcessor;

public class KitagawaProcessor implements PMCMCProcessor<KitagawaParams> {

	private List<String> logZs = new ArrayList<>();
	
	@Override
	public void output(File file) {
		PrintWriter writer = BriefIO.output(file);
		for (String line : logZs)
		{
			writer.println(line);
		}
		writer.close();
	}

	@Override
	public void process(KitagawaParams p) {
		logZs.add(p.asCommaSeparatedLine());
	}

}
