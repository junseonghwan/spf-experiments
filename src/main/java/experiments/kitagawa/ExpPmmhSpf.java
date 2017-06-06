package experiments.kitagawa;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import models.Kitagawa;

import org.apache.commons.lang3.tuple.Pair;

import pmcmc.LogZProcessor;
import pmcmc.MCMCProblemSpecification;
import pmcmc.PMCMCOptions;
import pmcmc.PMCMCProcessor;
import pmcmc.PMMHAlgorithm;
import simplesmc.AbstractSMCAlgorithm;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpPmmhSpf implements Runnable 
{
	@Option(required=false) public Random random = new Random(1);
	
	// parameters for data generation
	@Option(required=false) public double var_v = 10.0;
	@Option(required=false) public double var_w = 1.0;
	@Option(required=false) public int R = 500;
	
	// model hyperparameters
	@Option(required=false) public double shape = 1/0.01;
	@Option(required=false) public double rate = 1/0.01;
	@Option(required=false) public double sd_v = 0.15;
	@Option(required=false) public double sd_w = 0.08;
	
	// options
	@Option(name="pmcmc") public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@Option(name="spf") public SPFOptions spfOptions = new SPFOptions();

	@Override
	public void run() 
	{
		Pair<double[], double[]> ret = Kitagawa.simulate(random, var_v, var_w, R);

		List<PMCMCProcessor<KitagawaParams>> processors = new ArrayList<PMCMCProcessor<KitagawaParams>>();
		processors.add(new KitagawaProcessor("spf"));
		LogZProcessor<KitagawaParams> logZProcessor = new LogZProcessor<>("spf");

		MCMCProblemSpecification<KitagawaParams> mcmcProblemSpecification = new KitagawaMCMCProblemSpecification(shape, rate, sd_v, sd_w);
		KitagawaParams params = mcmcProblemSpecification.initialize(pmcmcOptions.random);
		spfOptions.numberOfConcreteParticles = 100;
		spfOptions.maxNumberOfVirtualParticles = 10000;
		spfOptions.targetedRelativeESS = 1.0;
		spfOptions.verbose = false;
		AbstractSMCAlgorithm<Double> spfAlgorithm = new StreamingParticleFilter<>(new KitagawaSMCProblemSpecification(params, ret.getRight()), spfOptions);
		PMMHAlgorithm<KitagawaParams, Double> pmmh = new PMMHAlgorithm<KitagawaParams, Double>(params, spfAlgorithm, mcmcProblemSpecification, pmcmcOptions, processors, logZProcessor);
		pmcmcOptions.nIter = 1000;
		pmcmcOptions.burnIn = 0;
		pmcmcOptions.thinningPeriod = 20;
		pmmh.sample(); // calling this function will generate the outputs
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpPmmhSpf());
	}

}
