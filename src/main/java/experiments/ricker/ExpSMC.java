package experiments.ricker;

import java.io.File;
import java.util.List;
import java.util.Random;

import models.RickerModel;

import org.apache.commons.lang3.tuple.Pair;

import simplesmc.GenericParticleProcessor;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;
import util.OutputHelper;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;

public class ExpSMC implements Runnable
{
	@Option(required=false) public static Random random = new Random(1721);
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.7;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 100;
	@Option(required=false) public static int numParticles = 10000;
	@Option(required=false) public static int numSimulations = 100;

	@Override
	public void run()
	{
		// generate the latent variable and the data (all at once)
		List<Pair<List<Double>, List<Integer>>> data = RickerModel.generate(random, numSimulations, T, N0, phi, var, r);

		// now, run SPF and infer the latent variables
		for (int i = 0; i < numSimulations; i++)
		{
			Pair<List<Double>, List<Integer>> ret = data.get(i);
			RickerParams params = new RickerParams(phi, r, var);

			SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, params, ret.getRight());
			SMCOptions options = new SMCOptions();
			options.nParticles = numParticles;
			options.random = new Random(random.nextLong());

			File outputDest = new File(Results.getResultFolder(), "/ricker-smc/simul" + (i+1) + "/");
			GenericParticleProcessor rickerParticleProcessor = new GenericParticleProcessor(outputDest.getAbsolutePath(), "particles/population");
			SMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(problemSpec, options, rickerParticleProcessor);
			smcAlgorithm.sample();

			OutputHelper.writeVector(new File(outputDest, "ricker-smc-latent.csv"), ret.getLeft());
			OutputHelper.writeVector(new File(outputDest, "ricker-smc-data.csv"), ret.getRight());
			OutputHelper.writeVector(new File(outputDest, "ricker-smc-ess.csv"), smcAlgorithm.effectiveSampleSize());
			OutputHelper.writeVector(new File(outputDest, "ricker-smc-timing.csv"), smcAlgorithm.timeInSeconds());

		}

		
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSMC());
	}
}
