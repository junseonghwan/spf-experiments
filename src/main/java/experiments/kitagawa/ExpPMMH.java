package experiments.kitagawa;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import models.Kitagawa;
import pmcmc.MCMCProblemSpecification;
import pmcmc.PMCMCOptions;
import pmcmc.PMCMCProcessor;
import pmcmc.PMMHAlgorithm;
import simplesmc.AbstractSMCAlgorithm;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpPMMH implements Runnable
{
	@Option(required=false) public Random random = new Random(1);
	@Option(required=false) public double var_v = 10.0;
	@Option(required=false) public double var_w = 1.0;
	@Option(required=false) public int R = 500;
	@Option(required=false) public double shape = 1/0.01;
	@Option(required=false) public double rate = 1/0.01;
	@Option(required=false) public double sd_v = 0.15;
	@Option(required=false) public double sd_w = 0.08;
	@Option(name="pmcmc") public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@Option(name="smc") public SMCOptions smcOptions = new SMCOptions();
	
	@Override
	public void run()
	{
		Pair<double[], double[]> ret = Kitagawa.simulate(random, var_v, var_w, R);
		
		List<PMCMCProcessor<KitagawaParams>> processors = new ArrayList<PMCMCProcessor<KitagawaParams>>();
		processors.add(new KitagawaProcessor());
		
		MCMCProblemSpecification<KitagawaParams> mcmcProblemSpecification = new KitagawaMCMCProblemSpecification(shape, rate, sd_v, sd_w);
		KitagawaParams params = mcmcProblemSpecification.initialize(pmcmcOptions.random);
		smcOptions.nParticles = 1000;
		AbstractSMCAlgorithm<Double> smcAlgorithm = new SMCAlgorithm<>(new KitagawaSMCProblemSpecification(params, ret.getRight()), smcOptions);
		PMMHAlgorithm<KitagawaParams, Double> pmmh = new PMMHAlgorithm<KitagawaParams, Double>(params, smcAlgorithm, mcmcProblemSpecification, pmcmcOptions, processors);
		pmcmcOptions.nIter = 1000;
		pmcmcOptions.burnIn = 0;
		pmcmcOptions.thinningPeriod = 20;
		pmmh.sample(); // calling this function will generate the outputs
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpPMMH());
	}

}
