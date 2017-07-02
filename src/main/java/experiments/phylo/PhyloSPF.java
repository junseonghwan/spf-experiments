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
import spf.SPFOptions;
import spf.StreamingParticleFilter;
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

public class PhyloSPF implements Runnable 
{
	@Option(required=false)
	public int numConcreteParticles = 10000;
	@Option(required=false)
	public int maxVirtualParticles = 100000;
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

		SPFOptions options = new SPFOptions();
		options.maxNumberOfVirtualParticles = maxVirtualParticles;
		options.numberOfConcreteParticles = numConcreteParticles;
		options.targetedRelativeESS = Double.POSITIVE_INFINITY;
		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		StreamingParticleFilter<PartialCoalescentState> spf = new StreamingParticleFilter<>(proposal, options);
		ParticlePopulation<PartialCoalescentState> pop = spf.sample();

		System.out.println(spf.logNormEstimate());

		// process the population
		double [][] dd = new double[numTaxa][numTaxa];
		double [][] trueDistances = new double[numTaxa][numTaxa];
		Counter<UnorderedPair<Taxon, Taxon>> truth = phylogeny.getPairwiseDistances(taxonIndexer);
		double trueHeight = phylogeny.getHeight();
		System.out.println(phylogeny.getTreeString());
		for (UnorderedPair<Taxon, Taxon> pair : truth)
		{
			int i = taxonIndexer.o2i(pair.getFirst());
			int j = taxonIndexer.o2i(pair.getSecond());
			double dist = truth.getCount(pair);
			trueDistances[i][j] = dist;
			trueDistances[j][i] = dist;
		}

		List<Double> heights = new ArrayList<>();
		for (PartialCoalescentState particle : pop.particles)
		{
			// compute the distance between two species, normalized by the height
			Counter<UnorderedPair<Taxon, Taxon>> distances = particle.getCoalescent().getPairwiseDistances(taxonIndexer);
			double h = particle.getCoalescent().getHeight();
			//System.out.println(h);
			heights.add(h);
			for (UnorderedPair<Taxon, Taxon> pair : distances)
			{
				int i = taxonIndexer.o2i(pair.getFirst());
				int j = taxonIndexer.o2i(pair.getSecond());
				double dist = distances.getCount(pair);
				dd[i][j] += dist;
				dd[j][i] += dist;
			}
		}
		for (int i = 0; i < numTaxa; i++)
		{
			for (int j = 0; j < numTaxa; j++)
			{
				dd[i][j] /= pop.particles.size();
			}
		}
		System.out.println("truth: " + phylogeny.getHeight());
		String [] header = new String[numTaxa];
		for (int i = 1; i <= numTaxa; i++) {
			header[i-1] = "T" + i;
		}
		OutputHelper.writeTableAsCSV(new File("output/phylo-pairwise-dist-spf.csv"), header, dd);
		OutputHelper.writeTableAsCSV(new File("output/phylo-pairwise-truth-spf.csv"), header, trueDistances);		
		OutputHelper.writeVector(new File("output/phylo-spf-ess.csv"), spf.relESS());
		OutputHelper.writeVector(new File("output/phylo-spf-heights.csv"), heights);
	}

	public static void main(String[] args) 
	{
		Mains.instrumentedRun(args, new PhyloSPF());
	}

}
