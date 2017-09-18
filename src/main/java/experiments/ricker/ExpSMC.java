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
import briefj.opt.OptionSet;
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
	@Option(required=false) public static int numSimulations = 1;
	//public final static String data_path = "output/ricker-data.csv";

	@OptionSet(name="smc")
	public SMCOptions option = new SMCOptions();
	
	@Override
	public void run()
	{
		// generate the latent variable and the data
		List<Pair<List<Double>, List<Integer>>> data = RickerModel.generate(random, numSimulations, T, N0, phi, var, r);
		// read the data
		//Pair<List<Double>, List<Integer>> ret = RickerModel.readFromFile(data_path);

		RickerModel model = new RickerModel(N0, phi, var, r);

		// now, run SMC and infer the latent variables
		for (int i = 0; i < numSimulations; i++)
		{
			Random rand = new Random(random.nextLong());
			option.random = rand;

			Pair<List<Double>, List<Integer>> ret = data.get(i);

			SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, model, ret.getRight());

			String outputDir = "/ricker-smc/";
			if (option.essThreshold == 0.0) 
				outputDir = "/ricker-sis/";

			File outputDest = new File(Results.getResultFolder(), outputDir + "simul" + (i+1) + "/");
			GenericParticleProcessor rickerParticleProcessor = new GenericParticleProcessor(outputDest.getAbsolutePath(), "population/");
			SMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(problemSpec, option, rickerParticleProcessor);
			smcAlgorithm.sample();

			OutputHelper.writeVector(new File(outputDest, "ricker-latent.csv"), ret.getLeft());
			OutputHelper.writeVector(new File(outputDest, "ricker-data.csv"), ret.getRight());
			OutputHelper.writeVector(new File(outputDest, "ricker-ess.csv"), smcAlgorithm.effectiveSampleSize());
			OutputHelper.writeVector(new File(outputDest, "ricker-timing.csv"), smcAlgorithm.timeInSeconds());
		}
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSMC());
	}
}
