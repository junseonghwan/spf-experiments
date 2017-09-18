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
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;

public class ExpPmmhSpf implements Runnable 
{
	@Option(required=false) public Random random = new Random(1);
	
	// parameters for data generation
	@Option(required=false) public double var_v = 10.0;
	@Option(required=false) public double var_w = 1.0;
	@Option(required=false) public int R = 500;
	
	// model hyperparameters
	@Option(required=false) public double sd_v = 0.15;
	@Option(required=false) public double sd_w = 0.08;
	
	// options
	@OptionSet(name="pmcmc") public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@OptionSet(name="spf") public SPFOptions spfOptions = new SPFOptions();

	@Override
	public void run() 
	{
		Pair<double[], double[]> ret = KitagawaModel.simulate(random, var_v, var_w, R);
		KitagawaModel model = new KitagawaModel(var_v, var_w);
		MultivariateUniformPrior prior = new MultivariateUniformPrior(new double[]{0.0, 10.0}, 2, false, true);

		LogZProcessor<RealVectorParameters> logZProcessor = new LogZProcessor<>("spf");

		MultivariateIndependentGaussianRandomWalk migrw = new MultivariateIndependentGaussianRandomWalk(new double[]{random.nextDouble(), random.nextDouble()}, new double[]{0.01, 0.01});
		AbstractSMCAlgorithm<Double> spfAlgorithm = new StreamingParticleFilter<>(new KitagawaSMCProblemSpecification(model, ret.getRight()), spfOptions);
		PMMHAlgorithm<RealVectorParameters, Double> pmmh = new PMMHAlgorithm<>(model, spfAlgorithm, migrw, prior, pmcmcOptions, null, logZProcessor, true);
		pmmh.sample(); // calling this function will generate the outputs
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpPmmhSpf());
	}

}
