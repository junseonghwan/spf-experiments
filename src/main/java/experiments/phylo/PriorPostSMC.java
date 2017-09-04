package experiments.phylo;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import bayonet.smc.ParticlePopulation;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;
import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.RootedPhylogenyProcessor;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantorModel;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import util.OutputHelper;

public class PriorPostSMC implements Runnable {

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
		
		for (int i = 0; i < numSimulations; i++)
			simulation(new Random(rand.nextLong()), i+1);
	}
	
	public double simulation(Random random, int simulNo)
	{

		// simulate a tree and generate the data
		EvolutionaryModel model = new JukesCantorModel(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);

		SMCOptions options = new SMCOptions();
		options.nParticles = numParticles;
		options.verbose = true;
		options.essThreshold = 1.0;
		options.random = new Random(random.nextLong());
		PriorPostProblemSpecification proposal = new PriorPostProblemSpecification(leaves);
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

		RootedPhylogenyProcessor processor = new RootedPhylogenyProcessor(phylogeny, taxonIndexer);
		int N = pop.particles.size();
		List<Double> heights = new ArrayList<>();
		System.out.println("Begin processing the particle population. Size=" + N);
		for (PartialCoalescentState particle : pop.particles)
		{
			processor.process(particle, 1.0/N);
			heights.add(particle.getCoalescent().getHeight());

		}
		Counter<UnorderedPair<Taxon, Taxon>> dist = processor.getMeanPairwiseDistances();
		for (UnorderedPair<Taxon, Taxon> key : dist.keySet())
		{
			int i = taxonIndexer.o2i(key.getFirst());
			int j = taxonIndexer.o2i(key.getSecond());
			dd[i][j] = dist.getCount(key);
			dd[j][i] = dd[i][j];
		}

		// compute the likelihood of the data given the true tree: p(y | t, \theta)
		System.out.println(phylogeny.getTreeString());
		FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, phylogeny);
		double logLik = PhyloOptions.calc.computeLoglik(phylogeny.getTaxon().getLikelihoodTable());
		System.out.println("p(y|t): " + logLik);
		System.out.println("truth: " + phylogeny.getHeight());
		// get the estimate of the marginal log likelihood
		double logZ = smc.logNormEstimate();
		System.out.println("logZ:" + logZ);

		String [] header = new String[numTaxa];
		for (int i = 1; i <= numTaxa; i++) {
			header[i-1] = "T" + i;
		}
		File resultsDir = Results.getResultFolder();
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-spf.csv"), header, dd);
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-spf.csv"), header, trueDistances);		
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-heights.csv"), heights);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-height-truth.csv"), new double[]{trueHeight});
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-logZ.csv"), new double[]{logZ});
		return smc.logNormEstimate();
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new PriorPostSMC());
	}

}