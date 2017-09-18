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
import briefj.opt.OptionSet;
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
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import util.OutputHelper;

public class TestPriorPostSMC implements Runnable {

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
			simulation(new Random(rand.nextLong()), i+1);
	}
	
	public double simulation(Random random, int simulNo)
	{

		// simulate a tree and generate the data
		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);

		options.random = new Random(random.nextLong());
		options.verbose = true;
		PriorPostProblemSpecification proposal = new PriorPostProblemSpecification(leaves);
		SMCAlgorithm<PartialCoalescentState> smc = new SMCAlgorithm<>(proposal, options);
		long start = System.currentTimeMillis();
		ParticlePopulation<PartialCoalescentState> pop = smc.sample();
		long end = System.currentTimeMillis();
		double samplingTime = (end - start)/1000.0;
		System.out.println(smc.logNormEstimate());

		// process the population
		if (phylogeny != null) {
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
			for (int j = 0; j < pop.particles.size(); j++)
			{
				PartialCoalescentState particle = pop.particles.get(j);
				processor.process(particle, pop.getNormalizedWeight(j));
				// record the heights
				heights.add(particle.getCoalescent().getHeight());
				weights.add(pop.getNormalizedWeight(j));

				// compute the log likelihood
				FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, particle.getCoalescent());
				double logLik = PhyloOptions.calc.computeLoglik(particle.getCoalescent().getTaxon().getLikelihoodTable());
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
			System.out.println("Processing time=" + particleProcessingTime);
	
			// compute the likelihood of the data given the true tree: p(y | t, \theta)
			System.out.println(phylogeny.getTreeString());
			FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, phylogeny);

			// TODO: compute the true log-likelihood via sum-product algorithm
			double logLik = PhyloOptions.calc.computeLoglik(phylogeny.getTaxon().getLikelihoodTable());
			System.out.println("p(y|t): " + logLik);
			System.out.println("truth: " + phylogeny.getHeight());
			// get the estimate of the marginal log likelihood
			double logZ = smc.logNormEstimate();
	
			String [] header = new String[numTaxa];
			for (int i = 0; i < numTaxa; i++) {
				header[i] = "T" + i;
			}
			//System.out.println(phylogeny.getNewickFormat());
			//System.out.println(phylogeny.getDataString());

			File resultsDir = Results.getResultFolder();
			OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-smc.csv"), header, dd);
			OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-smc.csv"), header, trueDistances);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-ess.csv"), smc.effectiveSampleSize());
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-heights.csv"), heights);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-weights.csv"), weights);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-height-truth.csv"), new double[]{trueHeight});
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-logZ.csv"), new double[]{logLik, logZ});
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-particle-logZs.csv"), logLiks);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-timing.csv"), new double[]{samplingTime, particleProcessingTime});
		} else {
			// just get the heights and the logZ
			List<Double> heights = new ArrayList<>();
			int N = pop.particles.size();
			System.out.println("Begin processing the particle population. Size=" + N);
			List<Double> logLiks = new ArrayList<>();
			for (PartialCoalescentState particle : pop.particles)
			{
				// record the heights
				heights.add(particle.getCoalescent().getHeight());
	
				// compute the log likelihood
				FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, particle.getCoalescent());
				double logLik = PhyloOptions.calc.computeLoglik(particle.getCoalescent().getTaxon().getLikelihoodTable());
				logLiks.add(logLik);
			}
			double logZ = smc.logNormEstimate();
			File resultsDir = Results.getResultFolder();
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-heights.csv"), heights);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-smc-logZ.csv"), new double[]{logZ});

		}
		return smc.logNormEstimate();

	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new TestPriorPostSMC());
	}

}
