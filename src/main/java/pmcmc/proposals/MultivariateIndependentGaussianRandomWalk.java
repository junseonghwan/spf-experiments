package pmcmc.proposals;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.MultivariateSummaryStatistics;

import bayonet.distributions.Normal;
import pmcmc.MCMCProposal;

public class MultivariateIndependentGaussianRandomWalk implements MCMCProposal<RealVectorParameters> {

	private double [] mu;
	private double [] sd;
	private int dim;
	
	public MultivariateIndependentGaussianRandomWalk(int dim, double [] sd) 
	{
		this.dim = dim;
		this.sd = sd;
		this.mu = new double[dim];
	}
	
	public MultivariateIndependentGaussianRandomWalk(double [] mu, double [] sd) 
	{
		this.dim = mu.length;
		this.sd = sd;
		this.mu = mu;
	}
	
	@Override
	public RealVectorParameters initialize(Random random) {
		return propose(random, new RealVectorParameters(mu));
	}

	@Override
	public RealVectorParameters propose(Random random, RealVectorParameters curr) {
		double [] vec = new double[dim];
		double [] currVec = curr.getVector();
		for (int i = 0; i < dim; i++) {
			vec[i] = Normal.generate(random, currVec[i], Math.pow(sd[i],2));
		}
		return new RealVectorParameters(vec);
	}

	@Override
	public double logProposalDensity(RealVectorParameters curr, RealVectorParameters prev) {
		double logProposalDensity = 0.0;
		double [] currVec = curr.getVector();
		double [] prevVec = prev.getVector();
		for (int i = 0; i < dim; i++) {
			logProposalDensity += Normal.logDensity(currVec[i], prevVec[i], Math.pow(sd[i],2));
		}
		return logProposalDensity;
	}

	@Override
	public void adapt(List<RealVectorParameters> params) {
		// compute the sample covariance matrix
		MultivariateSummaryStatistics summ = new MultivariateSummaryStatistics(dim, true);
		for (RealVectorParameters param : params)
		{
			summ.addValue(param.getVector());
		}
		sd = summ.getStandardDeviation();
		System.out.println("Adapting covariance matrix: ");
		for (int i = 0; i < dim; i++)
		{
			sd[i] *= Math.pow(2.38,2)/dim;
			System.out.println("sd[i]=" + sd[i]);
		}
	}
}
