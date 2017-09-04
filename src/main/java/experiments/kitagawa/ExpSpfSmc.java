package experiments.kitagawa;

import java.io.File;
import java.util.*;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.KitagawaModel;
import simplesmc.GenericParticleProcessor;
import simplesmc.ParticleProcessor;
import simplesmc.SMCProblemSpecification;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import spf.ResamplingScheme;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import util.OutputHelper;
import bayonet.distributions.Normal;
import bayonet.smc.ParticlePopulation;
import briefj.opt.Option;
import briefj.run.Mains;

public class ExpSpfSmc implements Runnable
{
	@Option(required=false) public static Random random = new Random(1);
	@Option(required=false) public static double var_v = 10.0;
	@Option(required=false) public static double var_w = 1.0;
	@Option(required=false) public static int R = 500;
	@Option(required=false) public static int numParticles = 1000;
	@Option(required=false) public static int numConcreteParticles = 1000;
	@Option(required=false) public static int maxVirtualParticles = 1000000;
	@Option(required=false) public static String smcOutputPath = "output/kitagawaSMC.csv";
	@Option(required=false) public static String spfOutputPath = "output/kitagawaSPF.csv";
	@Option(required=false) public static String generatedDataOutputPath = "output/kitagawaTruth.csv";
	@Option(required=false) public static boolean generateOutput = true;

	public static KitagawaSMCProblemSpecification proposal;

	@Override
	public void run()
	{
		experiment1();
	}

	/**
	 * compare SMC vs SPF with same number of concrete particles
	 */
	public void experiment1()
	{
		// 1. simulate from Kitagawa (1996)
		// 2. perform inference using SMC using K particles
		// 3. perform inference using SPF using K concrete particles
		// 4. compare the two distributions
		Pair<double [], double []> ret = KitagawaModel.simulate(random, var_v, var_w, R);
		double [] x = ret.getLeft();
		double [] y = ret.getRight();

		KitagawaModel model = new KitagawaModel(var_v, var_w);
		proposal = new KitagawaSMCProblemSpecification(model, y);
		
		OutputHelper.writeTableAsCSV(new File(generatedDataOutputPath), new String[]{"x", "y"}, x, y);

		// run SMC
		SMCOptions options = new SMCOptions();
		options.nParticles = numParticles;
		options.essThreshold = 1.0; // resample at every iteration
		ParticlePopulation<Double> population = runSMC(y, var_w, var_v, options);
		population = population.resample(random, bayonet.smc.ResamplingScheme.MULTINOMIAL);
		List<Double> particles = population.particles;
		// output to a file for plotting
		OutputHelper.writeVector(new File(smcOutputPath), particles);
		
		// run SPF
		List<Double> emissions = new ArrayList<>();
		for (double val : y) emissions.add(val);

		SPFOptions spfOptions = new SPFOptions();
		spfOptions.maxNumberOfVirtualParticles = maxVirtualParticles;
		spfOptions.numberOfConcreteParticles = numConcreteParticles;
		spfOptions.resamplingScheme = ResamplingScheme.MULTINOMIAL;
		spfOptions.targetedRelativeESS  = 1.0;
		spfOptions.verbose = false;
		ParticlePopulation<Double> samples = runSPF(emissions, var_w, var_v, spfOptions);
		OutputHelper.writeVector(new File(spfOutputPath), samples.particles);
	}
	
	private ParticlePopulation<Double> runSPF(List<Double> y, double var_w, double var_v, SPFOptions options)
	{		
		SMCProblemSpecification<Double> problemSpec = new SMCProblemSpecification<Double>() {

			@Override
			public Pair<Double, Double> proposeNext(int currentSmcIteration, Random random, Double currentParticle) {
				double mu = currentParticle/2.0 + 25*currentParticle / (1 + Math.pow(currentParticle, 2.0)) + 8*Math.cos(1.2*currentSmcIteration);
				double xstar = Normal.generate(random, mu, var_v);
				return Pair.of(logDensity(xstar, y.get(currentSmcIteration)), xstar);
		}

			@Override
			public Pair<Double, Double> proposeInitial(Random random) {
				double xstar = Normal.generate(random, 0.0, 5.0);
				double logw = Normal.logDensity(xstar, 0.0, 5.0);
				return Pair.of(logw, xstar);
			}

			@Override
			public int nIterations() {
				return y.size();
			}
			
			public double logDensity(Double latent, Double emission) {
				double logw = Normal.logDensity(emission, Math.pow(latent,  2.0)/20.0, var_w);
				return logw;
			}
		};
		
		ParticleProcessor<Double> particleProcessor = new GenericParticleProcessor("output/kitagawa-spf", "particles");
		StreamingParticleFilter<Double> spf = new StreamingParticleFilter<>(problemSpec, options, particleProcessor);
		ParticlePopulation<Double> population= spf.sample();
		return population;
	}

	private ParticlePopulation<Double> runSMC(double [] y, double var_w, double var_v, SMCOptions options)
	{
		ParticleProcessor<Double> particleProcessor = new GenericParticleProcessor("output/kitagawa-smc", "particles");
		SMCAlgorithm<Double> smc = new SMCAlgorithm<>(proposal, options, particleProcessor);
		return smc.sample(); // NOTE: resampling is not performed for the last iteration of SMC
	}
	
	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new ExpSpfSmc());
	}

}
