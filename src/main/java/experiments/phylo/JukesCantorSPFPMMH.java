package experiments.phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import briefj.Indexer;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningSystBiol2012;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantorModel;
import pmcmc.LogZProcessor;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import pmcmc.prior.MultivariateUniformPrior;
import pmcmc.proposals.MultivariateIndependentGaussianRandomWalk;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCProblemSpecification;
import spf.SPFOptions;
import spf.StreamingParticleFilter;

/**
 * Experiments to estimate the parameters of the evolutionary model.
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class JukesCantorSPFPMMH implements Runnable
{
	@Option(required=true)
	public Random rand = new Random(1);
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numTaxa = 4;
	@OptionSet(name="pmcmc")
	public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@OptionSet(name="spf")
	public SPFOptions spfOptions = new SPFOptions();
	@Option(required=true)
	public double trueMutationRate = 0.74;
	@Option(required=true)
	public double proposalVar = 0.1;
	@Option(required=true)
	public String proposalType = "priorprior";

	Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	List<Taxon> leaves = new ArrayList<Taxon>();

	private void generateData() {
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		// instantiate the model to generate the data and the tree
		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(trueMutationRate);
		PhyloOptions.calc = new FelsensteinPruningSystBiol2012(model);
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(rand, leaves);
		GenerateSequences.generateSequencesFromModel(rand, model, phylogeny, numSites);
	}
	
	@Override
	public void run()
	{
		generateData();

		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(rand.nextDouble()*5);
		PhyloOptions.calc = new FelsensteinPruningSystBiol2012(model);

		// pmcmc random
		pmcmcOptions.random = new Random(rand.nextLong());

		// specify SPF algorithm
		spfOptions.mainRandom = new Random(rand.nextLong());
		spfOptions.resamplingRandom = new Random(rand.nextLong());
		spfOptions.targetedRelativeESS = Double.POSITIVE_INFINITY;
		spfOptions.verbose = false;

		SMCProblemSpecification<PartialCoalescentState> proposal = null;
		if (proposalType.equalsIgnoreCase("priorprior"))
			proposal = new PriorPriorProblemSpecification(leaves);
		else if (proposalType.equalsIgnoreCase("priorpost"))
			proposal = new PriorPostProblemSpecification(leaves);

		LogZProcessor<RealVectorParameters> logZProcessor = new LogZProcessor<>("spf");

		StreamingParticleFilter<PartialCoalescentState> spf = new StreamingParticleFilter<>(proposal, spfOptions);

		MultivariateIndependentGaussianRandomWalk igrwProposal = new MultivariateIndependentGaussianRandomWalk(new double[]{2.5}, new double[]{proposalVar});
		MultivariateUniformPrior uniformPrior = new MultivariateUniformPrior(new double[]{0.0, 5.0}, 1, false, true);
		
		PMMHAlgorithm<RealVectorParameters, PartialCoalescentState> pmmh = new PMMHAlgorithm<>(model, spf, igrwProposal, uniformPrior, pmcmcOptions, null, logZProcessor, true);
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
		Mains.instrumentedRun(args, new JukesCantorSPFPMMH());
	}
}
