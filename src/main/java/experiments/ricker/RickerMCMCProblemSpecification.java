package experiments.ricker;

import java.util.Random;

import bayonet.distributions.Normal;
import bayonet.distributions.Uniform;
import pmcmc.MCMCProblemSpecification;

/**
 * Class representing a prior and random walk proposal on Ricker model parameters. 
 * The proposals are Gaussian distribution with variance 2.4, 1.2, 0.2
 * The prior is flat for each of the three parameters
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class RickerMCMCProblemSpecification implements MCMCProblemSpecification<RickerParams> 
{

	@Override
	public RickerParams initialize(Random random) {
		// randomly generate the parameters from Uniform prior
		double phi = Uniform.generate(random, 0, 20);
		double r = Math.exp(Uniform.generate(random, 0, 5));
		double var = Uniform.generate(random, 0.0, 1.0);
		RickerParams params = new RickerParams(phi, r, var);
		return params;
	}

	@Override
	public RickerParams propose(Random random, RickerParams curr) {
		// random walk proposal from Gaussian distribution with variance 1.5 for each of the parameters
		double phi = Normal.generate(random, curr.phi.getValue(), 1);
		double r = Normal.generate(random, curr.r.getValue(), 1);
		double var = Normal.generate(random, curr.var.getValue(), 1);
		RickerParams newParams = new RickerParams(phi, r, var);
		return newParams;
	}

	@Override
	public double logProposalDensity(RickerParams curr, RickerParams prev) {
		// because we are using the symmetric proposal, the logProposalDensity cancel
		return 0;
	}

	@Override
	public double logPriorDensity(RickerParams param) {
		if (param.phi.getValue() <= 0.0 || param.phi.getValue() >= 20.0)
			return Double.NEGATIVE_INFINITY;
		if (param.var.getValue() < 0.0)
			return Double.NEGATIVE_INFINITY;
		if (param.r.getValue() <= 0.0)
			return Double.NEGATIVE_INFINITY;
		return 0.0;
		//double logPrior = Exponential.logDensity(param.r.getValue(), 3.5);
		//logPrior += InverseGamma.logDensity(param.var.getValue(), 1.0, 1.0);
		//return logPrior;
	}

}
