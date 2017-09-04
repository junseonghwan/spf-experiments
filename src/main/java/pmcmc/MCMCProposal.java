package pmcmc;

import java.util.List;
import java.util.Random;

public interface MCMCProposal<P>
{
	public P initialize(Random random);
	public P propose(Random random, P curr);
	public double logProposalDensity(P curr, P prev);
	public void adapt(List<P> params);
}
