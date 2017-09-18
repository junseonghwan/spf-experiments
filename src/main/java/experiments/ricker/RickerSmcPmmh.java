package experiments.ricker;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.RickerModel;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import pmcmc.prior.MultivariateUniformPrior;
import pmcmc.proposals.MultivariateGaussianRandomWalk;
import pmcmc.proposals.MultivariateIndependentGaussianRandomWalk;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;
import util.OutputHelper;
import briefj.BriefFiles;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;

public class RickerSmcPmmh 
implements Runnable
{
	@Option(required=true) public static Random rand = new Random(1);
	@Option(required=false) public static double phi = 10.0;
	@Option(required=false) public static double r = 44.70118;
	@Option(required=false) public static double var = 0.3;
	@Option(required=false) public static double N0 = 7.0;
	@Option(required=false) public static int T = 50;
	@OptionSet(name="smc") public SMCOptions options = new SMCOptions();
	@OptionSet(name="pmcmc") public PMCMCOptions pmcmcOptions = new PMCMCOptions();

	@Override
	public void run()
	{
		Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(rand, T, N0, phi, var, r);

		double [] initial = new double[]{rand.nextDouble()*20, rand.nextDouble(), rand.nextDouble()*50};
		RickerModel model = new RickerModel(initial[0], initial[1], initial[2]);
		//RickerModel model = new RickerModel(phi, var, r);
		MultivariateUniformPrior prior = new MultivariateUniformPrior(new double[]{0.0, Double.POSITIVE_INFINITY}, 3, false, false);
		//double sqrtd = Math.sqrt(3.0);
		//double [] sd = new double[]{0.1/sqrtd, 0.1/sqrtd, 0.1/sqrtd};
		double [] sd = new double[]{1, 0.01, 1};
		//MultivariateIndependentGaussianRandomWalk mcmcProposal = new MultivariateIndependentGaussianRandomWalk(model.getModelParameters().getVector(), sd);
		MultivariateGaussianRandomWalk mcmcProposal = new MultivariateGaussianRandomWalk(model.getModelParameters().getVector(), sd);
		SMCProblemSpecification<Double> problemSpec = new RickerSMCProblemSpecification(T, model, ret.getRight());
		options.random = new Random(rand.nextLong());
		options.verbose = false;
		pmcmcOptions.random = new Random(rand.nextLong());
		SMCAlgorithm<Double> smc = new SMCAlgorithm<>(problemSpec, options);
		PMMHAlgorithm<RealVectorParameters, Double> pmmh = new PMMHAlgorithm<>(model, smc, mcmcProposal, prior, pmcmcOptions, true);
		pmmh.sample();
		//System.out.println(pmmh.nAccepts());
		File results = Results.getResultFolder();
		OutputHelper.writeVector(new File(results, "initialValues.csv"), initial);
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new RickerSmcPmmh());
	}

}
