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
import phylo.RootedPhylogenyProcessor;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantorModel;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import util.OutputHelper;
import bayonet.smc.ParticlePopulation;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;

public class TestPriorPriorSMC implements Runnable
{
	@Option(required=true)
	public Random rand = new Random(1172);
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numTaxa = 4;
	@Option(required=true)
	public int numSimulations = 50;
	@Option(required=true)
	public double mutationRate = 1.7;
	@OptionSet(name="smc")
	public SMCOptions options = new SMCOptions();

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
		{
			simulation(new Random(rand.nextLong()), i+1);			
		}
	}

	public double simulation(Random random, int simulNo)
	{
		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);

		options.random = new Random(random.nextLong());
		options.verbose = true;
		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		SMCAlgorithm<PartialCoalescentState> smc = new SMCAlgorithm<>(proposal, options);
		long start = System.currentTimeMillis();
		ParticlePopulation<PartialCoalescentState> pop = smc.sample();
		long end = System.currentTimeMillis();
		double samplingTime = (end - start)/1000.0;
		double logZ = smc.logNormEstimate();

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
		List<Double> weights = new ArrayList<>();
		System.out.println("Begin processing the particle population. Size=" + N);
		start = System.currentTimeMillis();
		List<Double> logLiks = new ArrayList<>();

		for (int l = 0; l < pop.particles.size(); l++)
		{
			PartialCoalescentState state = pop.particles.get(l);
			double ww = pop.getNormalizedWeight(l);
			processor.process(state, ww);
			weights.add(ww);
			heights.add(state.getCoalescent().getHeight());

			// compute the log likelihood
			FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, state.getCoalescent());
			double logLik = PhyloOptions.calc.computeLoglik(state.getCoalescent().getTaxon().getLikelihoodTable());
			logLiks.add(logLik);
		}
		Counter<UnorderedPair<Taxon, Taxon>> dist = processor.getMeanPairwiseDistances();
		for (UnorderedPair<Taxon, Taxon> key : dist.keySet())
		{
			int i = taxonIndexer.o2i(key.getFirst());
			int j = taxonIndexer.o2i(key.getSecond());
			dd[i][j] = dist.getCount(key);
			dd[j][i] = dd[i][j];
		}

		end = System.currentTimeMillis();
		double particleProcessingTime = (end - start)/1000.0;

		// compute the likelihood of the data given the true tree: p(y | t, \theta)
		System.out.println(phylogeny.getTreeString());
		FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, phylogeny);
		double dataLogLik = PhyloOptions.calc.computeLoglik(phylogeny.getTaxon().getLikelihoodTable());
		System.out.println("p(y|t)=" + dataLogLik);
		double logZR = smc.getLogNorms().get(leaves.size() - 2);

		System.out.println("truth: " + phylogeny.getHeight());
		String [] header = new String[numTaxa];
		for (int i = 1; i <= numTaxa; i++) 
		{
			header[i-1] = "T" + i;
		}
		File resultsDir = Results.getResultFolder();
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-smc.csv"), header, dd);
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-smc.csv"), header, trueDistances);		
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-ess.csv"), smc.effectiveSampleSize());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-heights.csv"), heights);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-weights.csv"), weights);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-height-truth.csv"), new double[]{trueHeight});
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-logZ.csv"), new double[]{dataLogLik, logZR, logZ});
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-particle-logZs.csv"), logLiks);
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-timing.csv"), new double[]{samplingTime, particleProcessingTime});
		return smc.logNormEstimate();
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new TestPriorPriorSMC());
	}
}
