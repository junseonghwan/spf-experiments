package experiments.ricker;

import java.util.List;
import java.util.Random;

import models.RickerModel;

import org.apache.commons.lang3.tuple.Pair;

import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpSMC implements Runnable
{
	@Option(required=false) public static Random random = new Random(721);
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.7;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 10;
	@Option(required=false) public static int numParticles = 1000;

	@Override
	public void run()
	{
		Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(random, T, N0, phi, var, r);
		RickerParams params = new RickerParams(phi, r, var);

		SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, params, ret.getRight());
		SMCOptions options = new SMCOptions();
		options.nParticles = numParticles;
		SMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(problemSpec, options);
		smcAlgorithm.sample();
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSMC());
	}
}
