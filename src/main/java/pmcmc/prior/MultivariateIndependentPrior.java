package pmcmc.prior;

import java.util.List;

import bayonet.distributions.UnivariateRealDistribution;
import pmcmc.ProbabilityDistribution;
import pmcmc.proposals.RealVectorParameters;

public class MultivariateIndependentPrior implements ProbabilityDistribution<RealVectorParameters> 
{
	private List<Class<UnivariateRealDistribution>> distributions;
	public MultivariateIndependentPrior(List<Class<UnivariateRealDistribution>> distributions)
	{
		this.distributions = distributions;
	}

	@Override
	public double logDensity(RealVectorParameters p) {
		double [] vec = p.getVector();
		double logDensity = 0.0;
		for (int i = 0; i < vec.length; i++)
		{
			
		}
		return 0;
	}

}
