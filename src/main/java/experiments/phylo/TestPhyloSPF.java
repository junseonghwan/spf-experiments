package experiments.phylo;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.RootedPhylogenyProcessor;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantor;
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import util.OutputHelper;
import bayonet.smc.ParticlePopulation;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;

public class TestPhyloSPF implements Runnable 
{
	@Option(required=true)
	public double targetESS = 1.0;
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
	public int numSimulations = 20;

	public static Set<Integer> simulationSet = new HashSet<>();
	private Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	private List<Taxon> leaves = new ArrayList<Taxon>();

	//public static String data_path = "output/DNA_data.txt";
	public static String data_path = null;
	/*
    public static CTMCExpFam<String> model = CTMCExpFam.createModelWithFullSupport(Indexers.dnaIndexer(), true);
    static {
	    model.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<String>()));
	    model.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<String>()));
    }
    */
	
	static {
		int [] simulationNo = new int[]{};
		for (int i : simulationNo)
		{
			simulationSet.add(i);
		}
	}

	@Override
	public void run() 
	{
		if (data_path != null) {
			int n = 0;
			for (String line : BriefIO.readLines(new File(data_path)))
			{
				Taxon T = new Taxon("T" + n++);
				T.setSequence(line.trim());
				leaves.add(T);
				taxonIndexer.addToIndex(T);
			}
			simulation(new Random(rand.nextLong()), 1, false);
		} else {
			for (int n = 0; n < numTaxa; n++)
			{
				Taxon T = new Taxon("T" + n);
				leaves.add(T);
				taxonIndexer.addToIndex(T);
			}


			for (int i = 0; i < numSimulations; i++)
			{
				System.out.println("Simulation: " + (i+1));
				Random random = new Random(rand.nextLong());
				if (simulationSet.size() == 0)
					simulation(random, i+1, true);
				else if (simulationSet.size() > 0 && simulationSet.contains(i+1))
					simulation(random, i+1, true);
			}
		}
	}
	
	private double simulation(Random random, int simulNo, boolean generateData)
	{
		// simulate a tree and generate the data
		/*
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

		EvolutionaryModel model = new JukesCantor(1.0);
		/*
		double [] pi = new double[]{0.3,0.2,0.2,0.3};
		//GTRModelParams gtrParams = new GTRModelParams(rand.nextDouble(), rand.nextDouble(), rand.nextDouble(), rand.nextDouble(), rand.nextDouble(), rand.nextDouble());
		GTRModelParams gtrParams = new GTRModelParams(0.26,0.18,0.17,0.15,0.11,0.13);
		EvolutionaryModel model = new GTRModel(pi, gtrParams);
		System.out.println(gtrParams.toString());
		*/

		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		RootedPhylogeny phylogeny = null;
		if (generateData) {
			phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
			//GenerateSequenceDataFromExpFam.generateSequencesCTMC(random, model, learnedModel, phylogeny, numSites);
			GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);
		}

		SPFOptions options = new SPFOptions();
		options.maxNumberOfVirtualParticles = maxVirtualParticles;
		options.numberOfConcreteParticles = numConcreteParticles;
		options.mainRandom = new Random(random.nextLong());
		options.resamplingRandom = new Random(random.nextLong());
		options.targetedRelativeESS = Double.POSITIVE_INFINITY;
		options.storeParticleWeights = true;
		//options.targetedRelativeESS = targetESS;
		options.verbose = true;
		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		StreamingParticleFilter<PartialCoalescentState> spf = new StreamingParticleFilter<>(proposal, options);
		long start = System.currentTimeMillis();
		ParticlePopulation<PartialCoalescentState> pop = spf.sample();
		long end = System.currentTimeMillis();
		double samplingTime = (end - start)/1000.0;
		System.out.println("Sampling time=" + samplingTime);

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
			System.out.println("Begin processing the particle population. Size=" + N);
			start = System.currentTimeMillis();
			List<Double> logLiks = new ArrayList<>();
			for (PartialCoalescentState particle : pop.particles)
			{
				processor.process(particle, 1.0/N);
				// record the heights
				heights.add(particle.getCoalescent().getHeight());
	
				// compute the log likelihood
				FelsensteinPruningAlgorithm.computeDataLogLikTable((FelsensteinPruningAlgorithm)PhyloOptions.calc, particle.getCoalescent());
				double logLik = PhyloOptions.calc.computeLoglik(particle.getCoalescent().getTaxon().getLikelihoodTable());
				logLiks.add(logLik);
				
				// compute the distance between two species
				/*
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
				*/
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
			double logZ = spf.logNormEstimate();
			double logZR = spf.getLogZs().get(leaves.size() - 2);
			System.out.println("logZ:" + logZ + ", logZR:" + logZR);
	
			String [] header = new String[numTaxa];
			for (int i = 0; i < numTaxa; i++) {
				header[i] = "T" + i;
			}
			System.out.println(phylogeny.getNewickFormat());
			System.out.println(phylogeny.getDataString());

			File resultsDir = Results.getResultFolder();
			OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-spf.csv"), header, dd);
			OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/phylo-pairwise-dist-truth-spf.csv"), header, trueDistances);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-ess.csv"), spf.relESS());
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-heights.csv"), heights);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-height-truth.csv"), new double[]{trueHeight});
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-logZ.csv"), new double[]{logLik, logZR, logZ});
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-particle-logZs.csv"), logLiks);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-timing.csv"), new double[]{samplingTime, particleProcessingTime});
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
			double logZ = spf.logNormEstimate();
			double logZR = spf.getLogZs().get(leaves.size() - 2);
			File resultsDir = Results.getResultFolder();
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-heights.csv"), heights);
			OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/phylo-spf-logZ.csv"), new double[]{logZR, logZ});

		}
		return spf.logNormEstimate();
	}

	public static void main(String[] args) 
	{
		Mains.instrumentedRun(args, new TestPhyloSPF());
	}

}
