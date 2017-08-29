package experiments.phylo;

import java.util.Random;

import bayonet.distributions.Normal;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PhyloOptions;
import phylo.models.JukesCantor;
import phylo.models.JukesCantor.JukesCantorParam;
import pmcmc.MCMCProblemSpecification;

public class JukesCantorMCMCProblemSpec implements MCMCProblemSpecification<JukesCantorParam> 
{
	private JukesCantorParam curr;
	private double var;

	public JukesCantorMCMCProblemSpec(double var) {
		this.var = var;
	}
	
	@Override
	public JukesCantorParam initialize(Random random) 
	{
		if (curr == null)
			curr = new JukesCantorParam(random.nextDouble());
		return curr;
	}

	@Override
	public JukesCantorParam propose(Random random, JukesCantorParam curr) 
	{
		if (curr.getMutationRate() <= 0.0)
			throw new RuntimeException();
		// do a random walk 
		double delta = Normal.generate(random, 0.0, var);
		double newValue = delta + curr.getMutationRate();
		JukesCantorParam newParam = new JukesCantorParam(newValue);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(new JukesCantor(newValue));
		return newParam;
	}

	@Override
	public double logProposalDensity(JukesCantorParam curr, JukesCantorParam prev) 
	{
		if (curr.getMutationRate() <= 0.0)
			return Double.NEGATIVE_INFINITY;
		
		// using Normal distribution so the ratio cancels
		return 0;
	}

	@Override
	public double logPriorDensity(JukesCantorParam param) 
	{
		// let's assume Uniform prior on the parameter from (0,infty)
		if (param.getMutationRate() <= 0.0)
			return Double.NEGATIVE_INFINITY;
		return 0;
	}

}
