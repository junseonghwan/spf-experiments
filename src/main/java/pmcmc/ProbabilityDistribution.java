package pmcmc;

public interface ProbabilityDistribution<P extends ModelParameters> 
{
	public double logDensity(P p);
}
