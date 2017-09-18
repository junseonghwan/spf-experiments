package simplesmc;

import java.io.File;

import bayonet.smc.ParticlePopulation;
import util.OutputHelper;

public class GenericParticleProcessor implements ParticleProcessor<Double> 
{
	private String outputFilePrefix;
	private String outputLocation;
	
	public GenericParticleProcessor(String outputLocation, String outputFilePrefix) 
	{
		this.outputLocation = outputLocation;
		this.outputFilePrefix = outputFilePrefix;
	}

	@Override
	public void process(int currentIteration, ParticlePopulation<Double> population) 
	{
		//File results = Results.getResultFolder();
		double [] weights = new double[population.nParticles()];
		for (int i = 0; i < population.nParticles(); i++)
			weights[i] = population.getNormalizedWeight(i);
		OutputHelper.writeVector(new File(outputLocation + "/" + outputFilePrefix + "weights" + currentIteration + ".csv"), weights);
		OutputHelper.writeVector(new File(outputLocation + "/" + outputFilePrefix + "particles" + currentIteration + ".csv"), population.particles);
	}

}
