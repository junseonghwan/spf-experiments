package experiments.phylo;

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
import spf.SPFOptions;
import spf.StreamingParticleFilter;
import briefj.Indexer;
import briefj.opt.Option;
import briefj.run.Mains;

public class JukesCantorPriorPrioSPF implements Runnable {

	@Option(required=true)
	public Random rand = new Random(1);
	@Option(required=true)
	public int numConcreteParticles = 100;
	@Option(required=true)
	public int maxVirtualParticles = 10000;
	@Option(required=true)
	public double targetESS = 1.0;
	@Option(required=true)
	public int numTaxa = 20;
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numSimulations = 5;
	@Option(required=true)
	public double mutationRate = 1.0;

	private Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	private List<Taxon> leaves = new ArrayList<Taxon>();
	private RootedPhylogeny phylogeny = null;

	@Override
	public void run() 
	{
		// generate the data once, run simulations multiple times	
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		EvolutionaryModel model = new JukesCantor(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		Random random = new Random(rand.nextLong());
		phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);

		for (int i = 0; i < numSimulations; i++)
		{
			random = new Random(rand.nextLong());
			System.out.println("Simulation: " + (i+1));
			simulation(random, i+1);
		}
	}
	
	public void simulation(Random random, int simulNo)
	{
		SPFOptions options = new SPFOptions();
		options.maxNumberOfVirtualParticles = maxVirtualParticles;
		options.numberOfConcreteParticles = numConcreteParticles;
		options.mainRandom = new Random(random.nextLong());
		options.resamplingRandom = new Random(random.nextLong());
		options.storeParticleWeights = true;
		options.targetedRelativeESS = targetESS;
		options.verbose = true;

		PriorPriorProblemSpecification proposal = new PriorPriorProblemSpecification(leaves);
		StreamingParticleFilter<PartialCoalescentState> spf = new StreamingParticleFilter<>(proposal, options);
		spf.sample();
		PartialCoalescentStateProcessorUtil.generateOutputSPF(simulNo, taxonIndexer, spf, phylogeny);
	}

	public static void main(String[] args) 
	{
		Mains.instrumentedRun(args, new JukesCantorPriorPrioSPF());
	}

}
