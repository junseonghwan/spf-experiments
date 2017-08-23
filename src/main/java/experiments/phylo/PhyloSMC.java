package experiments.phylo;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantor;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import util.OutputHelper;
import bayonet.smc.ParticlePopulation;
import briefj.BriefParallel;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;

public class PhyloSMC implements Runnable
{
	@Option(required=true)
	public int numParticles = 1000;
	@Option(required=true)
	public Random rand = new Random(1172);
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numTaxa = 4;
	@Option(required=true)
	public int numSimulations = 100;
	@Option(required=true)
	public double mutationRate = 0.7;

	private Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	private List<Taxon> leaves = new ArrayList<Taxon>();

	@Override
	public void run()
	{
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		BriefParallel.process(numSimulations, 4, i -> {simulation(new Random(rand.nextLong()), i+1);});
	}

	public double simulation(Random random, int simulNo)
	{
		// simulate a tree and generate the data
		/*
	    CTMCExpFam<String> model = CTMCExpFam.createModelWithFullSupport(Indexers.dnaIndexer(), true);
	    model.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<String>()));
	    model.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<String>()));

		int p = model.nFeatures();

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
		LikelihoodCalculatorExpFam calc = new LikelihoodCalculatorExpFam(model, learnedModel);
		*/

		EvolutionaryModel model = new JukesCantor(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		//GenerateSequenceDataFromExpFam.generateSequencesCTMC(random, model, learnedModel, phylogeny, numSites);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);

		SMCOptions options = new SMCOptions();
		options.nParticles = numParticles;
		options.verbose = false;
		options.random = new Random(random.nextLong());
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
			/*
			System.out.println(particle.toString());
			System.out.println(ww);
			*/
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
		
		// compute the likelihood of the data given the true tree: p(y | t, \theta)
		System.out.println(phylogeny.getTreeString());
		FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, phylogeny);
		System.out.println(PhyloOptions.calc.computeLoglik(phylogeny.getTaxon().getLikelihoodTable()));

		System.out.println("truth: " + phylogeny.getHeight());
		String [] header = new String[numTaxa];
		for (int i = 1; i <= numTaxa; i++) {
			header[i-1] = "T" + i;
		}
		File resultsDir = Results.getResultFolder();
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-smc.csv"), header, dd);
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-smc.csv"), header, trueDistances);		
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-ess.csv"), smc.effectiveSampleSize());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-heights.csv"), heights);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-weights.csv"), weights);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-height-truth.csv"), new double[]{trueHeight});

		return smc.logNormEstimate();
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new PhyloSMC());
	}
}
