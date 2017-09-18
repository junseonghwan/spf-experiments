package experiments.phylo;

import briefj.run.Mains;
import phylo.PhyloUtils;

public class PrimatesDataAnalysis implements Runnable 
{

	@Override
	public void run() 
	{
		PhyloUtils.readData("data/primates.txt");
	}

	public static void main(String[] args) 
	{
		Mains.instrumentedRun(args, new PrimatesDataAnalysis());
	}

}
