package pmcmc;

import java.util.Random;

public interface MCMCProblemSpecification<P> 
{
	public P initialize(Random random);
	public P propose(Random random, P curr);
	public double logProposalDensity(P curr, P prev);
	public double logPriorDensity(P param);
}
