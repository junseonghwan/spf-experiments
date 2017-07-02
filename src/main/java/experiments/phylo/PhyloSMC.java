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
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;

public class PhyloSMC implements Runnable
{
	@Option(required=false)
	public int numParticles = 10000;
	@Option(required=false)
	public Random rand = new Random(172);
	@Option(required=false)
	public int numSites = 1000;
	@Option(required=false)
	public int numTaxa = 10;
	@Option(required=false)
	public double rate = 1.2;
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
		PhyloOptions.rate = rate;

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
		
		// process the particles
		double [][] dd = new double[numTaxa][numTaxa];
		double [][] trueDistances = new double[numTaxa][numTaxa];
		Counter<UnorderedPair<Taxon, Taxon>> truth = phylogeny.getPairwiseDistances(taxonIndexer);
		double trueHeight = phylogeny.getHeight();
		for (UnorderedPair<Taxon, Taxon> pair : truth)
		{
			int i = taxonIndexer.o2i(pair.getFirst());
			int j = taxonIndexer.o2i(pair.getSecond());
			double dist = truth.getCount(pair);
			trueDistances[i][j] = dist;
			trueDistances[j][i] = dist;
		}

		List<Double> heights = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		for (int l = 0; l < pop.particles.size(); l++)
		{
			PartialCoalescentState particle = pop.particles.get(l);
			double ww = pop.getNormalizedWeight(l);
			System.out.println(particle.toString());
			System.out.println(ww);
			weights.add(ww);
			// compute the distance between two species, normalized by the height
			Counter<UnorderedPair<Taxon, Taxon>> distances = particle.getCoalescent().getPairwiseDistances(taxonIndexer);
			double h = particle.getCoalescent().getHeight();
			heights.add(h);
			for (UnorderedPair<Taxon, Taxon> pair : distances)
			{
				int i = taxonIndexer.o2i(pair.getFirst());
				int j = taxonIndexer.o2i(pair.getSecond());
				double dist = distances.getCount(pair);
				dd[i][j] += ww*dist;
				dd[j][i] += ww*dist;
			}
		}
		/*
		for (int i = 0; i < numTaxa; i++)
		{
			for (int j = 0; j < numTaxa; j++)
			{
				dd[i][j] /= pop.particles.size();
			}
		}
		*/
		
		// compute the likelihood of the truth
		System.out.println(phylogeny.getTreeString());
		System.out.println(PhyloOptions.calc.computeLogLik(phylogeny));
		
		System.out.println("truth: " + phylogeny.getHeight());
		String [] header = new String[numTaxa];
		for (int i = 1; i <= numTaxa; i++) {
			header[i-1] = "T" + i;
		}
		OutputHelper.writeTableAsCSV(new File("output/phylo-pairwise-dist-smc.csv"), header, dd);
		OutputHelper.writeTableAsCSV(new File("output/phylo-pairwise-truth-smc.csv"), header, trueDistances);		
		OutputHelper.writeVector(new File("output/phylo-smc-ess.csv"), smc.effectiveSampleSize());
		OutputHelper.writeVector(new File("output/phylo-smc-heights.csv"), heights);
		OutputHelper.writeVector(new File("output/phylo-smc-weights.csv"), weights);
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new PhyloSMC());
	}
}
