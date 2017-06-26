package experiments.ricker;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import simplesmc.SMCProblemSpecification;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import models.RickerModel;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpSPF implements Runnable
{
	@Option(required=false) public static Random random = new Random(721);
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.7;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 10;
	@Option(required=false) public static int numConcreteParticles = 1000;
	@Option(required=false) public static int maxVirtualParticles = 10000;

	@Override
	public void run()
	{
		Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(random, T, N0, phi, var, r);
		RickerParams params = new RickerParams(phi, r, var);

		SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, params, ret.getRight());
		SPFOptions options = new SPFOptions();
		options.maxNumberOfVirtualParticles = maxVirtualParticles;
		options.numberOfConcreteParticles = numConcreteParticles;
		options.targetedRelativeESS = 1.0;
		StreamingParticleFilter<Double> spf = new StreamingParticleFilter<>(problemSpec, options);
		spf.sample();
		//ParticlePopulation<Double> population= spf.sample();
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSPF());
	}

}