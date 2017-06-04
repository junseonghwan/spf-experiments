package simplesmc;

import bayonet.smc.ParticlePopulation;

public abstract class AbstractSMCAlgorithm<P> 
{
	public abstract ParticlePopulation<P> sample();
	public abstract double logNormEstimate();
}
