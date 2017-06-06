package simplesmc;

import java.util.List;

import bayonet.smc.ParticlePopulation;

public abstract class AbstractSMCAlgorithm<P> 
{
	protected List<Double> timeInSeconds;
	
	public abstract ParticlePopulation<P> sample();
	public abstract double logNormEstimate();
	public List<Double> timeInSeconds() { return timeInSeconds; }
}
