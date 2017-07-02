package experiments.ricker;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import simplesmc.GenericParticleProcessor;
import simplesmc.SMCProblemSpecification;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import util.OutputHelper;
import models.RickerModel;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;

public class ExpSPF implements Runnable
{
	@Option(required=false) public static Random random = new Random(1721);
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.7;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 100;
	@Option(required=false) public static int numConcreteParticles = 1000;
	@Option(required=false) public static int maxVirtualParticles = 10000;
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
			SPFOptions options = new SPFOptions();
			options.maxNumberOfVirtualParticles = maxVirtualParticles;
			options.numberOfConcreteParticles = numConcreteParticles;
			options.targetedRelativeESS = 1.0;
			options.resamplingRandom = new Random(random.nextLong());
			options.mainRandom = new Random(random.nextLong());

			File outputDest = new File(Results.getResultFolder(), "/ricker-spf/simul" + (i+1) + "/");
			GenericParticleProcessor rickerParticleProcessor = new GenericParticleProcessor(outputDest.getAbsolutePath(), "particles/population");
			StreamingParticleFilter<Double> spf = new StreamingParticleFilter<>(problemSpec, options, rickerParticleProcessor);
			spf.sample();

			OutputHelper.writeVector(new File(outputDest, "ricker-spf-latent.csv"), ret.getLeft());
			OutputHelper.writeVector(new File(outputDest, "ricker-spf-data.csv"), ret.getRight());
			OutputHelper.writeVector(new File(outputDest, "ricker-spf-nimplicit.csv"), spf.nImplicitParticles());
			OutputHelper.writeVector(new File(outputDest, "ricker-spf-reless.csv"), spf.relESS());
			OutputHelper.writeVector(new File(outputDest, "ricker-spf-timing.csv"), spf.timeInSeconds());
		}
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSPF());
	}

}
