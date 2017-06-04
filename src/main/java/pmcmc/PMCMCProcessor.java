package pmcmc;

import java.io.File;

public interface PMCMCProcessor<P> {
	public void process(P p);
	public void output(File file);
}
