package experiments.ricker;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.RickerModel;
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
	@Option(required=false) public static double r = 44.70118;
	@Option(required=false) public static double c = 1.0;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 50;
	@Option(required=false) public static int numParticles = 1000;
	@Option(required=false) public static int numSimulations = 1;
	public final static String data_path = "output/ricker-data.csv";

	@Override
	public void run()
	{
		// generate the latent variable and the data
		//List<Pair<List<Double>, List<Integer>>> data = RickerModel.generate(random, numSimulations, T, N0, phi, var, r);
		// read the data
		Pair<List<Double>, List<Integer>> ret = RickerModel.readFromFile(data_path);
		
		RickerModel model = new RickerModel(N0, phi, var, r);

		// now, run SMC and infer the latent variables
		for (int i = 0; i < numSimulations; i++)
		{
			Random rand = new Random(random.nextLong());
			//Pair<List<Double>, List<Integer>> ret = data.get(i);

			SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, model, ret.getRight());
			SMCOptions options = new SMCOptions();
			options.nParticles = numParticles;
			options.random = rand;

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
