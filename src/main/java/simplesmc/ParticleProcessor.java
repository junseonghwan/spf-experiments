package simplesmc;

import bayonet.smc.ParticlePopulation;

public interface ParticleProcessor<T> 
{
	public void process(int currentIteration, ParticlePopulation<T> population);
}
