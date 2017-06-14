package experiments.kitagawa;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import models.Kitagawa;
import pmcmc.LogZProcessor;
import pmcmc.MCMCProblemSpecification;
import pmcmc.PMCMCOptions;
import pmcmc.PMCMCProcessor;
import pmcmc.PMMHAlgorithm;
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
		Pair<double[], double[]> ret = Kitagawa.simulate(random, var_v, var_w, R);

		List<PMCMCProcessor<KitagawaParams>> processors = new ArrayList<PMCMCProcessor<KitagawaParams>>();
		processors.add(new KitagawaProcessor("smc"));
		LogZProcessor<KitagawaParams> logZProcessor = new LogZProcessor<>("smc");

		MCMCProblemSpecification<KitagawaParams> mcmcProblemSpecification = new KitagawaMCMCProblemSpecification(shape, rate, sd_v, sd_w);
		KitagawaParams params = mcmcProblemSpecification.initialize(pmcmcOptions.random);
		AbstractSMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(new KitagawaSMCProblemSpecification(params, ret.getRight()), smcOptions);
		PMMHAlgorithm<KitagawaParams, Double> pmmh = new PMMHAlgorithm<KitagawaParams, Double>(params, smcAlgorithm, mcmcProblemSpecification, pmcmcOptions, processors, logZProcessor);
		pmmh.sample(); // calling this function generates the output containing the parameters, logZ estimates, and number of acceptances
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpPmmhSmc());
	}

}
