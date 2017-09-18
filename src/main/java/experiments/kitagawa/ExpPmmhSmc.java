package experiments.kitagawa;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.KitagawaModel;
import pmcmc.LogZProcessor;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import pmcmc.prior.MultivariateUniformPrior;
import pmcmc.proposals.MultivariateIndependentGaussianRandomWalk;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.AbstractSMCAlgorithm;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;

public class ExpPmmhSmc implements Runnable
{
	@Option(required=false) public Random random = new Random(1);
	
	// parameters for data generation
	@Option(required=false) public double var_v = 10.0;
	@Option(required=false) public double var_w = 1.0;
	@Option(required=false) public int R = 500;

	// model hyperparameters
	@Option(required=false) public double shape = 0.01;
	@Option(required=false) public double rate = 0.01;
	@Option(required=false) public double sd_v = 0.15;
	@Option(required=false) public double sd_w = 0.08;
	
	// options
	@OptionSet(name="pmcmc") public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@OptionSet(name="smc") public SMCOptions smcOptions = new SMCOptions();
	
	@Override
	public void run()
	{
		Pair<double[], double[]> ret = KitagawaModel.simulate(random, var_v, var_w, R);

		LogZProcessor<RealVectorParameters> logZProcessor = new LogZProcessor<>("smc");

		KitagawaModel model = new KitagawaModel(var_v, var_w);
		MultivariateUniformPrior prior = new MultivariateUniformPrior(new double[]{0.0, 10.0}, 2, false, true);
		MultivariateIndependentGaussianRandomWalk migrw = new MultivariateIndependentGaussianRandomWalk(new double[]{random.nextDouble(), random.nextDouble()}, new double[]{0.01, 0.01});
		AbstractSMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(new KitagawaSMCProblemSpecification(model, ret.getRight()), smcOptions);
		PMMHAlgorithm<RealVectorParameters, Double> pmmh = new PMMHAlgorithm<>(model, smcAlgorithm, migrw, prior, pmcmcOptions, null, logZProcessor, true);
		pmmh.sample(); // calling this function generates the output containing the parameters, logZ estimates, and number of acceptances
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpPmmhSmc());
	}

}
