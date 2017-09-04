package pmcmc.prior;

import pmcmc.ProbabilityDistribution;
import pmcmc.proposals.RealVectorParameters;

public class MultivariateUniformPrior implements ProbabilityDistribution<RealVectorParameters> 
{
	private double [][] support;
	private boolean [] minSupportInclusive;
	private boolean [] maxSupportInclusive;
	public MultivariateUniformPrior(double [][] support, boolean [] minSupportInclusive, boolean [] maxSupportInclusive)
	{
		this.support = support;
		this.minSupportInclusive = minSupportInclusive;
		this.maxSupportInclusive = maxSupportInclusive;
	}
	public MultivariateUniformPrior(double [][] support, boolean minSupportInclusive, boolean maxSupportInclusive)
	{
		this.support = support;
		this.minSupportInclusive = new boolean[support.length];
		this.maxSupportInclusive = new boolean[support.length];
		for (int i = 0; i < support.length; i++) {
			this.maxSupportInclusive[i] = minSupportInclusive;
			this.maxSupportInclusive[i] = maxSupportInclusive;
		}
	}
	public MultivariateUniformPrior(double [] support, boolean minSupportInclusive, boolean maxSupportInclusive)
	{
		this.support = new double[support.length][2];
		this.minSupportInclusive = new boolean[support.length];
		this.maxSupportInclusive = new boolean[support.length];
		for (int i = 0; i < support.length; i++) {
			this.support[i] = support;
			this.maxSupportInclusive[i] = minSupportInclusive;
			this.maxSupportInclusive[i] = maxSupportInclusive;
		}
	}

	@Override
	public double logDensity(RealVectorParameters p) 
	{
		double [] vec = p.getVector();
		for (int i = 0; i < vec.length; i++) {
			if ((minSupportInclusive[i] && vec[i] < support[i][0]) ||
					(!minSupportInclusive[i] && vec[i] <= support[i][0]) ||
					(!maxSupportInclusive[i] && vec[i] > support[i][1]) ||
					(maxSupportInclusive[i] && vec[i] >= support[i][1])) 
				return Double.NEGATIVE_INFINITY;
		}
		return 0.0;
	}

}
