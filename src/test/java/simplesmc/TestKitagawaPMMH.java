package simplesmc;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import distributions.InverseGamma;
import experiments.kitagawa.KitagawaParams;
import experiments.kitagawa.KitagawaSMCProblemSpecification;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.variables.RealVariable;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import models.Kitagawa;
import simplesmc.SMCOptions;
import simplesmc.pmcmc.PMCMCFactor;
import spf.SPFOptions;
import spf.StreamingParticleFilter;

public class TestKitagawaPMMH implements Runnable {

	@Option(required=false) public static Random random = new Random(1);
	@Option(required=false) public static double var_v = 10.0;
	@Option(required=false) public static double var_w = 1.0;
	@Option(required=false) public static double shape = 0.01;
	@Option(required=false) public static double rate = 0.01;
	@Option(required=false) public static int R = 10;
	@OptionSet(name="smc") public static SMCOptions smcOptions = new SMCOptions();
	public static SPFOptions spfOptions = new SPFOptions();
	@OptionSet(name = "factory") public final MCMCFactory factory = new MCMCFactory();

	public static double [] observations;

	/**
	 * Model to be used for sampling the variance of the Nonlinear SSM of Kitagawa (1996)
	 * 
	 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
	 *
	 */
	public class Model
	{
		@FactorArgument(makeStochastic=false)
		public RealVariable a = new RealVariable(shape);
		@FactorArgument(makeStochastic=false)
		public RealVariable b = new RealVariable(rate);
		@FactorComponent
		public KitagawaParams params = new KitagawaParams();

	    @DefineFactor
	    //public PMCMCFactor<Double> likelihood = new PMCMCFactor<>(params, new SMCAlgorithm<>(new KitagawaProblemSpecification(params, observations), smcOptions));
	    public PMCMCFactor<Double> likelihood = new PMCMCFactor<>(params, new StreamingParticleFilter<>(new KitagawaSMCProblemSpecification(params, observations), spfOptions));
	    @DefineFactor
	    public InverseGamma<InverseGamma.ShapeRateParameterization> var_v = InverseGamma.on(params.var_v, a, b);
	    @DefineFactor
	    public InverseGamma<InverseGamma.ShapeRateParameterization> var_w = InverseGamma.on(params.var_w, a, b);
	}
	
	@Override
	public void run() 
	{
		// 1. generate the data
		// 2. sample the parameters
		
		Pair<double [], double []> ret = Kitagawa.simulate(random, var_v, var_w, R);
		observations = ret.getRight();
		
	    Model model = new Model();
	    MCMCAlgorithm mcmc = factory.build(model, false);
	    mcmc.options.nMCMCSweeps = 1000;
	    mcmc.options.burnIn = 100;
	    System.out.println(mcmc.model);
	    mcmc.run();
		
	}

	public static void main(String[] args) {
		Mains.instrumentedRun(args, new TestKitagawaPMMH());
	}

}
