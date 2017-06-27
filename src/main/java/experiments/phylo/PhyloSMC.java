package experiments.phylo;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import phylo.CustomRandomGenerator;
import phylo.GenerateData;
import phylo.LikelihoodCalculator;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import util.OutputHelper;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCExpFam.LearnedReversibleModel;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import conifer.io.Indexers;
import bayonet.smc.ParticlePopulation;
import briefj.Indexer;
import briefj.opt.Option;
import briefj.run.Mains;

public class PhyloSMC implements Runnable
{
	@Option(required=false)
	public int numParticles = 10000;
	@Option(required=false)
	public Random rand = new Random(1);
	@Option(required=false)
	public double rate = 1.71;
	@Option(required=false)
	public int numSites = 100;
	@Option(required=false)
	public int numTaxa = 4;
	@Option(required=false)
	public int numSimulations = 1;

	@Override
	public void run()
	{
		// simulate a tree and generate the data
	    CTMCExpFam<String> model = CTMCExpFam.createModelWithFullSupport(Indexers.dnaIndexer(), true);
	    model.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<String>()));
	    model.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<String>()));

		int p = model.nFeatures();
		System.out.println("nFeatures: " + p);

		// generate w ~ mvn(0, I)
		double [] mu = new double[p];
		double [][] I = new double[p][p];
		for (int i = 0; i < p; i++)
		{
			I[i][i] = 1.0;
		}

		MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new CustomRandomGenerator(rand.nextLong()), mu, I);
		double [] w = mvn.sample();

		@SuppressWarnings("rawtypes")
		LearnedReversibleModel learnedModel = model.reversibleModelWithParameters(w);
		LikelihoodCalculator calc = new LikelihoodCalculator(model, learnedModel);
		
		Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
		List<Taxon> leaves = new ArrayList<Taxon>();
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		PhyloOptions.calc = calc;

		// generate the data and the tree
		RootedPhylogeny phylogeny = GenerateData.sampleRootedPhylogeny(rand, leaves, PhyloOptions.rate);
		GenerateData.generateSequencesCTMC(rand, model, learnedModel, phylogeny, numSites);

		SMCOptions options = new SMCOptions();
		options.nParticles = numParticles;
		options.verbose = true;
		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		SMCAlgorithm<PartialCoalescentState> smc = new SMCAlgorithm<>(proposal, options);
		ParticlePopulation<PartialCoalescentState> pop = smc.sample();
		System.out.println(smc.logNormEstimate());
		
		// process the population
		List<Double> heights = new ArrayList<Double>();
		for (PartialCoalescentState particle : pop.particles)
		{
			heights.add(particle.getCoalescent().getHeight());
		}
		System.out.println("truth: " + phylogeny.getHeight());
		OutputHelper.writeVector(new File("output/phylo-smc.csv"), heights);
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new PhyloSMC());
	}
}
