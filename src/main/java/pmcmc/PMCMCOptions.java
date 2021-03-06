package pmcmc;

import java.util.Random;

import briefj.opt.Option;

public class PMCMCOptions 
{
	  @Option(gloss = "Seed for the MCMC algorithm")
	  public Random random = new Random(1);
	  
	  @Option(gloss = "Number of iterations")
	  public int nIter = 10000;

	  @Option(gloss = "Number of burnin interations")
	  public int burnIn = 100;

	  @Option(gloss = "Thinning period. Should be greater or equal to 1.") 
	  public int thinningPeriod = 20;
	  
	  public boolean verbose = false;

}
