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

		OutputHelper.writeVector(new File(outputLocation + "/" + outputFilePrefix + currentIteration + ".csv"), population.particles);
	}

}
