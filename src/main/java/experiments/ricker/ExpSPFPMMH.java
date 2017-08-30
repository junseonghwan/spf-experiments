package experiments.ricker;

import java.util.List;
import java.util.Random;

import models.RickerModel;

import org.apache.commons.lang3.tuple.Pair;

import pmcmc.MCMCProblemSpecification;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import simplesmc.SMCProblemSpecification;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpSPFPMMH 
implements Runnable
{
	@Option(required=false) public static Random random = new Random(721);
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.70118;
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 50;
	@Option(required=false) public static int numConcreteParticles = 100;
	@Option(required=false) public static int maxVirtualParticles = 1000;

	@Override
	public void run()
	{
		Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(random, T, N0, phi, var, r);
		MCMCProblemSpecification<RickerParams> mcmcProblemSpec = new RickerMCMCProblemSpecification();
		
		RickerParams params = new RickerParams(15, 20, 1);

		SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, params, ret.getRight());
		SPFOptions spfOptions = new SPFOptions();
		spfOptions.maxNumberOfVirtualParticles = maxVirtualParticles;
		spfOptions.numberOfConcreteParticles = numConcreteParticles;
		spfOptions.targetedRelativeESS = Double.POSITIVE_INFINITY;
		spfOptions.verbose = false;
		StreamingParticleFilter<Double> spf = new StreamingParticleFilter<>(problemSpec, spfOptions);
		
		PMCMCOptions pmcmcOptions = new PMCMCOptions();
		pmcmcOptions.burnIn = 1000;
		pmcmcOptions.nIter = 10000;
		pmcmcOptions.thinningPeriod = 1;
		PMMHAlgorithm<RickerParams, Double> pmmh = new PMMHAlgorithm<RickerParams, Double>(params, spf, mcmcProblemSpec, pmcmcOptions);
		pmmh.sample();
		System.out.println(pmmh.nAccepts());
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSPFPMMH());
	}

}
