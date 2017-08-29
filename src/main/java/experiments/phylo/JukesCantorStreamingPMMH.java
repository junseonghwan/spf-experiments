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
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GTRModel;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantor;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import phylo.models.GTRModel.GTRModelParams;
import phylo.models.JukesCantor.JukesCantorParam;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import util.OutputHelper;

/**
 * Experiments to estimate the parameters of the evolutionary model.
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class JukesCantorStreamingPMMH implements Runnable
{
	@Option(required=true)
	public int numConcreteParticles = 100;
	@Option(required=true)
	public int maxVirtualParticles = 10000;
	@Option(required=true)
	public Random rand = new Random(1);
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numTaxa = 4;
	@Option(required=true)
	public int chainLength = 10000;
	@Option(required=true)
	public double trueMutationRate = 0.74;
	@Option(required=true)
	public double var = 0.01;

	@Override
	public void run()
	{
		Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
		List<Taxon> leaves = new ArrayList<Taxon>();
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		// instantiate the model to generate the data and the tree
		EvolutionaryModel model = new JukesCantor(trueMutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(rand, leaves);
		GenerateSequences.generateSequencesFromModel(rand, model, phylogeny, numSites);

		// specify SPF algorithm
		SPFOptions options = new SPFOptions();
		options.maxNumberOfVirtualParticles = maxVirtualParticles;
		options.numberOfConcreteParticles = numConcreteParticles;
		options.mainRandom = new Random(rand.nextLong());
		options.resamplingRandom = new Random(rand.nextLong());
		options.targetedRelativeESS = Double.POSITIVE_INFINITY;
		options.verbose = true;
		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		StreamingParticleFilter<PartialCoalescentState> spf = new StreamingParticleFilter<>(proposal, options);

		JukesCantorParam params = new JukesCantorParam(1.2);
		model = new JukesCantor(params.getMutationRate());
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);
		JukesCantorMCMCProblemSpec mcmcProblemSpecification = new JukesCantorMCMCProblemSpec(var);
		PMCMCOptions pmcmcOptions = new PMCMCOptions();
		pmcmcOptions.random = new Random(rand.nextLong());
		pmcmcOptions.burnIn = 0;
		pmcmcOptions.nIter = chainLength;
		PMMHAlgorithm<JukesCantorParam, PartialCoalescentState> pmmh = new PMMHAlgorithm<>(params, spf, mcmcProblemSpecification, pmcmcOptions);
		pmmh.sample();
	}

//	private double simulation(Random random, int simulNo)
//	{
//
//		// process the population
//		double [][] dd = new double[numTaxa][numTaxa];
//		double [][] trueDistances = new double[numTaxa][numTaxa];
//		Counter<UnorderedPair<Taxon, Taxon>> truth = phylogeny.getPairwiseDistances(taxonIndexer);
//		double trueHeight = phylogeny.getHeight();
//		for (UnorderedPair<Taxon, Taxon> pair : truth)
//		{
//			int i = taxonIndexer.o2i(pair.getFirst());
//			int j = taxonIndexer.o2i(pair.getSecond());
//			double dist = truth.getCount(pair);
//			trueDistances[i][j] = dist;
//			trueDistances[j][i] = dist;
//		}
//
//		List<Double> heights = new ArrayList<>();
//		for (PartialCoalescentState particle : pop.particles)
//		{
//			// compute the distance between two species
//			Counter<UnorderedPair<Taxon, Taxon>> distances = particle.getCoalescent().getPairwiseDistances(taxonIndexer);
//			double h = particle.getCoalescent().getHeight();
//			//System.out.println(h);
//			heights.add(h);
//			for (UnorderedPair<Taxon, Taxon> pair : distances)
//			{
//				int i = taxonIndexer.o2i(pair.getFirst());
//				int j = taxonIndexer.o2i(pair.getSecond());
//				double dist = distances.getCount(pair);
//				dd[i][j] += dist;
//				dd[j][i] += dist;
//			}
//		}
//		for (int i = 0; i < numTaxa; i++)
//		{
//			for (int j = 0; j < numTaxa; j++)
//			{
//				dd[i][j] /= pop.particles.size();
//			}
//		}
//
//		// compute the likelihood of the data given the true tree: p(y | t, \theta)
//		System.out.println(phylogeny.getTreeString());
//		FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, phylogeny);
//		double logLik = PhyloOptions.calc.computeLoglik(phylogeny.getTaxon().getLikelihoodTable());
//		System.out.println("p(y|t): " + logLik);
//		System.out.println("truth: " + phylogeny.getHeight());
//		// get the estimate of the marginal log likelihood
//		double logZ = spf.logNormEstimate();
//		System.out.println("logZ:" + logZ);
//		
//		String [] header = new String[numTaxa];
//		for (int i = 1; i <= numTaxa; i++) {
//			header[i-1] = "T" + i;
//		}
//		File resultsDir = Results.getResultFolder();
//		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-spf.csv"), header, dd);
//		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-spf.csv"), header, trueDistances);		
//		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-ess.csv"), spf.relESS());
//		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-heights.csv"), heights);
//		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-height-truth.csv"), new double[]{trueHeight});
//		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-logZ.csv"), new double[]{numConcreteParticles, maxVirtualParticles, (end - start)/1000.0, logZ});
//		
//		return spf.logNormEstimate();
//	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new JukesCantorStreamingPMMH());
	}
}
